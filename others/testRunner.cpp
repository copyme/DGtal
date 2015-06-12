#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <vector>
#include "DGtal/curves/ParametricCurveDigitizer3D.h"
#include "DGtal/curves/parametric/EllipticHelix.h"
#include <DGtal/geometry/curves/SaturatedSegmentation.h>
#include "DGtal/curves/parametric/DecoratorCurveTransformation.h"
#include "DGtal/images/RigidTransformation3D.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/Naive3DDSSComputer.h"
#include "DGtal/geometry/curves/LambdaMST2D.h"
#include "DGtal/geometry/curves/LambdaMST3D.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include <DGtal/io/boards/Board3D.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/base/SimpleConstRange.h>
#include <fenv.h>
#include <climits>

using namespace std;
using namespace DGtal;
using namespace Z3i;
using namespace functors;

template <typename TCurve>
void print3D ( TCurve & digitalCurve,  Board3D<> & board )
{
  for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
  {
    if( Z3i::Point(97, 99, 278) != digitalCurve.at(i) )
    {
      board.setFillColor ( Color ( 0, 0, 255, 100 ) );
      board << SetMode3D(digitalCurve[i].className(), "");
      board << digitalCurve[i];
    }
    else
    {
      board.setFillColor ( Color ( 255, 0, 0, 100 ) );
      board << SetMode3D(digitalCurve[i].className(), "");
      board << digitalCurve[i];
      board << SetMode3D(digitalCurve[i].className(), "PavingWired");
      board << digitalCurve[i];
    }
  }
}

template <typename TCurve>
void print2D ( TCurve const & digitalCurve,  Board2D & board )
{
  for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
  {
    board.setFillColor ( Color ( 0, 0, 255, 100 ) );
    board << SetMode(digitalCurve[i].className(), "");
    board << digitalCurve[i];
    
    board << SetMode(digitalCurve[i].className(), "Both");
    board << digitalCurve[i];
  }
}

template <class TPoint>
bool areConnected ( TPoint const & a, TPoint const & b )
{
  for ( unsigned short i = 0; i < TPoint::dimension; i++ )
  {
    if ( std::abs ( a[i] - b[i] ) > 1 )
      return false;
  }
  return true;
}

template <class TPoint>
double distance ( TPoint const & a, TPoint const & b )
{
  double dist = 0.;
  for ( unsigned short i = 0; i < TPoint::dimension; i++ )
  {
    dist += ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  return std::sqrt ( dist );
}

template <class Iterator, class TPoint>
Iterator closest_point ( TPoint const & p, Iterator begin, Iterator end )
{
  double d = 0;
  Iterator result = begin;
  while (begin != end)
  {
    double d2 = distance( p, *begin );
    if (d2 < d)
    {
      d = d2;
      result = begin;
    }
    ++begin;
  }
  return result;
}

template <class TCurve, class TPoint >
void geometricalSort ( TCurve & curve )
{ 
  TCurve tmp;
  std::copy( curve.begin(), curve.end(), std::back_inserter ( tmp ) );
  curve.clear();
  
  TPoint p = tmp.front();
  tmp.erase( tmp.begin() );
  curve.push_back ( p );
  while (tmp.size())
    {
    // Find next point (the one closest to p)
    typename TCurve::iterator it = closest_point<typename TCurve::iterator, TPoint>( p, tmp.begin(), tmp.end() );
    curve.push_back( *it );
    // Our new p is the point we just found
    p = *it;
    // Remove the point we just found from the vector of points
    tmp.erase( it );
    }
}

template <typename Projector >
void project ( const Projector & proj, std::vector < Z3i::Point > & curve3d, std::vector < Z2i::Point > & curve2d )
{
  for ( typename std::vector<Z3i::Point>::const_iterator it = curve3d.begin(); it != curve3d.end(); ++it )
  {
    if ( curve2d.size() == 0 )
      curve2d.push_back ( proj ( *it ) );
    else if ( std::find ( curve2d.begin(), curve2d.end(), proj ( *it ) ) == curve2d.end() )
      curve2d.push_back ( proj ( *it ) );
  }
  geometricalSort< std::vector < Z2i::Point >, Z2i::Point > ( curve2d );
}

template <typename TCurve>
void print2DXY ( TCurve & digitalCurve,  Board3D<> & board, int distance )
{
  for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
  {
    board.setFillColor ( Color ( 255, 0, 0, 100 ) );
    board << SetMode3D(digitalCurve[i].className(), "");
    board << Z3i::Point ( digitalCurve.at(i)[0], digitalCurve.at(i)[1], distance );
    
    board << SetMode3D(digitalCurve[i].className(), "PavingWired");
    board << Z3i::Point ( digitalCurve.at(i)[0], digitalCurve.at(i)[1], distance );
  }
}

template <typename TCurve>
void print2DXZ ( TCurve & digitalCurve,  Board3D<> & board, int distance )
{
  for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
  {
    board.setFillColor ( Color ( 0, 255, 0, 100 ) );
    board << SetMode3D(digitalCurve[i].className(), "");
    board << Z3i::Point ( digitalCurve.at(i)[0], distance, digitalCurve.at(i)[1] );
    
    board << SetMode3D(digitalCurve[i].className(), "PavingWired");
    board << Z3i::Point ( digitalCurve.at(i)[0], distance, digitalCurve.at(i)[1] );
  }
}

template <typename TCurve>
void print2DYZ ( TCurve & digitalCurve,  Board3D<> & board, int distance )
{
  for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
  {
    board.setFillColor ( Color ( 0, 0, 255, 100 ) );
    board << SetMode3D(digitalCurve[i].className(), "");
    board << Z3i::Point ( distance, digitalCurve.at(i)[0], digitalCurve.at(i)[1] );
    
    board << SetMode3D(digitalCurve[i].className(), "PavingWired");
    board << Z3i::Point ( distance, digitalCurve.at(i)[0], digitalCurve.at(i)[1] );
  }
}

template <typename TSegmentComp>
void printDSS3D ( TSegmentComp & segmenter, Board3D<> & board )
{
  for ( typename TSegmentComp::SegmentComputerIterator it = segmenter.begin() ; it != segmenter.end(); ++it )
  {
    if ( (*it).isInDSS ( Z3i::Point(97, 99, 278) ) )
    {
      board << SetMode3D((*it).className(), "BoundingBox");
      board.setFillColor ( Color ( 255, 255, 255, 255 ) );
      board << *it;
    }
  }
}

template <typename TSegmentComp>
void printDSS2D ( TSegmentComp & segmenter, Board2D & board )
{
  for ( typename TSegmentComp::SegmentComputerIterator it = segmenter.begin() ; it != segmenter.end(); ++it )
  {
    board << SetMode("ArithmeticalDSS", "BoundingBox");
    board.setFillColor ( Color ( 255, 255, 255, 255 ) );
    board << it->primitive();
  }
}

template <typename TDCurve>
void exportPINK ( double radius1, double radius2, double depth, const TDCurve & digitalCurve )
{
  std::ofstream curve ( "curve" + std::to_string ( radius1 ) + "_" + std::to_string ( radius2 ) + "_" + std::to_string ( depth ) + ".list" ); // pink like format
  curve << "B " << digitalCurve.size() << std::endl;
  for (int i = 0; i  < digitalCurve.size(); i++ )
    curve << digitalCurve.at(i)[0] << " " << digitalCurve.at(i)[1] << " " << digitalCurve.at(i)[2] << std::endl;
  curve.close();
}

template <typename TCurve>
unsigned int findMainAxis ( const TCurve & curve, double t )
{
  Z3i::RealVector value;
  value[0] = std::fabs ( curve.xp ( t )[0] );
  value[1] = std::fabs ( curve.xp ( t )[1] );
  value[2] = std::fabs ( curve.xp ( t )[2] );
  
  if ( std::isgreater ( value[0], value[1] ) &&  std::isgreater ( value[0], value[2] ) )
    return 0;
  else if ( std::isgreater( value[1], value[0] ) && std::isgreater ( value[1], value[2] ) )
    return 1;
  else if ( std::isgreater ( value[2], value[0] ) && std::isgreater ( value[2], value[1] ) )
    return 2;
}

template <typename TCurve>
typename TCurve::const_iterator FindClosestPoint ( Z3i::RealPoint & x, TCurve & dCurve )
{
  typename TCurve::const_iterator closestIt = dCurve.begin();
  double minDist = DBL_MAX;
  for (typename TCurve::const_iterator it = dCurve.begin(); it != dCurve.end(); ++it )
  {
    Z3i::Point y = *it;
    double distance = std::sqrt ( ( (x[0] - y[0]) * (x[0] - y[0]) ) + ( (x[1] - y[1]) * (x[1] - y[1]) ) +  ( (x[2] - y[2]) * (x[2] - y[2]) ) );
    if ( minDist > distance )
    {
      minDist = distance;
      closestIt = it;
    }
  }
  return closestIt;
}

template <class TCurve, class TPoint>
bool isCloseCurve ( TCurve & curve )
{
  TPoint & start = curve.front();
  TPoint & end = curve.back();
  
  for ( unsigned short i = 0; i < TPoint::dimension; i++ )
  {
    if ( std::abs ( start[i] - end[i] ) > 1 )
      return false;
  }
  return true;
}

//split to smaller classes
template < class TCurve >
class TestRunner
{
  private: //types
    typedef ForwardRigidTransformation3D < Space, Identity, Space::RealPoint, Space::RealPoint > ForwardTrans;
    typedef BackwardRigidTransformation3D < Space, Identity, Space::RealPoint, Space::RealPoint > BackwardTrans;
    typedef DecoratorCurveTransformation < TCurve, ForwardTrans, BackwardTrans > MyRotatedCurve;
    typedef ParametricCurveDigitizer3DDecorator < MyRotatedCurve >  Digitizer;
    typedef typename ParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef std::vector< Z2i::Point > MyProjection;
    typedef std::vector < Z2i::Point > Range;
    typedef Range::const_iterator ConstIterator;
    typedef ArithmeticalDSSComputer < typename MyProjection::const_iterator, int, 8 > SegmentComputer;
    typedef ArithmeticalDSSComputer < SimpleConstRange<typename MyProjection::const_iterator>::ConstCirculator, int, 8 > SegmentComputerCloseCurve;
    typedef Naive3DDSSComputer < typename MyDigitalCurve::const_iterator, int, 8 > Naive3SegmentComputer;
    typedef Naive3DDSSComputer < typename SimpleConstRange<typename MyDigitalCurve::const_iterator>::ConstCirculator, int, 8 > Naive3SegmentComputerCloseCurve;
    typedef SaturatedSegmentation<SegmentComputer> Segmentation;
    typedef SaturatedSegmentation<SegmentComputerCloseCurve> SegmentationCloseCurve;
    typedef SaturatedSegmentation<Naive3SegmentComputer> SegmentationNaive3D;
    typedef SaturatedSegmentation<Naive3SegmentComputerCloseCurve> SegmentationNaive3DCloseCurve;
    typedef Projector< SpaceND<2,int> > Projector2d;
    
    private: //members
      float radius1, radius2, depth;
      double angle;
      Z3i::RealVector axis;
      double mse2d, mse3d, mae2d, mae3d;
      bool lengthFilter;
      bool distanceFilter;
      TCurve * realCurve;
      ForwardTrans * trans;
      BackwardTrans * inverse;
      MyRotatedCurve * rotCurve;
      MyDigitalCurve digitalCurve;
      SegmentationCloseCurve * segmenterXYCloseCurve;
      Segmentation * segmenterXY;
      SegmentationCloseCurve * segmenterXZCloseCurve;
      Segmentation * segmenterXZ;
      SegmentationCloseCurve * segmenterYZCloseCurve;
      Segmentation * segmenterYZ;
      Naive3SegmentComputerCloseCurve * computer3DCloseCurve;
      SegmentationNaive3DCloseCurve * segmenter3DCloseCurve;
      Naive3SegmentComputer * computer3D;
      SegmentationNaive3D * segmenter3D;
      LambdaMST3D < SegmentationNaive3DCloseCurve > * lmst3D64CloseCurve;
      LambdaMST3D < SegmentationNaive3D > * lmst3D64;
      LambdaMST2D < SegmentationCloseCurve > * lmst64CloseCurve;
      LambdaMST2D < Segmentation > * lmst64;
      MyProjection oXY;
      MyProjection oXZ;
      MyProjection oYZ;
      Projector2d myProjXY, myProjXZ, myProjYZ;
      unsigned int lenMin;
      Board3D<> & board;
      Board2D board2d;
      
      double mint, maxt;
      
      void isInitialized()
      {
	if ( realCurve == nullptr || trans == nullptr || inverse == nullptr || rotCurve == nullptr  )
	  throw "Not initialized! First call; setCurveParams; setTransformationParams.";
      }
      
      void generateProjections()
      {
	isInitialized();
	oXY.clear();
	oXZ.clear();
	oYZ.clear();
	project ( myProjXY, digitalCurve, oXY );
// 	print2D<MyProjection>( oXY, board2d );
	project ( myProjXZ, digitalCurve, oXZ );
	project ( myProjYZ, digitalCurve, oYZ );
	
// 	print2DXY< MyProjection >( oXY, board, -( std::max<double>(radius1, radius2) + std::min<double>( radius1, radius2)/2. ) );
// 	print2DXZ< MyProjection >( oXZ, board, -( std::max<double>(radius1, radius2) + std::min<double>( radius1, radius2)/2. ) );
// 	print2DYZ< MyProjection >( oYZ, board, -( std::max<double>(radius1, radius2) + std::min<double>( radius1, radius2)/2. ) );
      }
      
      void remove2DSegements()
      {
	delete segmenterXYCloseCurve;
	segmenterXYCloseCurve = nullptr;
	delete segmenterXY;
	segmenterXY = nullptr;
	delete segmenterXZCloseCurve;
	segmenterXZCloseCurve = nullptr;
	delete segmenterXZ;
	segmenterXZ = nullptr;
	delete segmenterYZCloseCurve;
	segmenterYZCloseCurve = nullptr;
	delete segmenterYZ;
	segmenterYZ = nullptr;
      }
      
      void remove3DSegements()
      {
	delete computer3DCloseCurve;
	computer3DCloseCurve = nullptr;
	delete segmenter3DCloseCurve;
	segmenter3DCloseCurve = nullptr;
	delete computer3D;
	computer3D = nullptr;
	delete segmenter3D;
	segmenter3D = nullptr;
      }
      
      void removeEstimators()
      {
	delete lmst3D64CloseCurve,
	lmst3D64CloseCurve = nullptr,
	delete lmst3D64;
	lmst3D64 = nullptr;
	delete lmst64CloseCurve;
	lmst64CloseCurve = nullptr;
	delete lmst64;
	lmst64 = nullptr;
      }
      
      void initEstimators()
      {
	remove3DSegements();
	removeEstimators();
	
	computer3D = new Naive3SegmentComputer();
	segmenter3D = new SegmentationNaive3D ( digitalCurve.begin(), digitalCurve.end(), *computer3D );
	
	lmst3D64CloseCurve = new LambdaMST3D < SegmentationNaive3DCloseCurve >();
	lmst3D64 = new LambdaMST3D < SegmentationNaive3D >();
	lmst64CloseCurve = new LambdaMST2D < SegmentationCloseCurve >();
	lmst64 = new  LambdaMST2D < Segmentation >();
	
	
	lmst3D64->init ( digitalCurve.begin(), digitalCurve.end() );
	lmst3D64->attach ( *segmenter3D );
	
	
	if ( isCloseCurve<MyDigitalCurve, Z3i::Point> ( digitalCurve ) )
	{
	  SimpleConstRange<typename MyDigitalCurve::const_iterator> t_Range(digitalCurve.begin(), digitalCurve.end());
	  computer3DCloseCurve = new Naive3SegmentComputerCloseCurve();
	  segmenter3DCloseCurve = new SegmentationNaive3DCloseCurve ( t_Range.c(), --t_Range.c(), *computer3DCloseCurve );
	  repairCover<MyDigitalCurve, SegmentationNaive3DCloseCurve>( digitalCurve, *segmenter3DCloseCurve );
	  lmst3D64CloseCurve->init ( t_Range.c(), --t_Range.c() );
	  lmst3D64CloseCurve->attach ( *segmenter3DCloseCurve );
	  printDSS3D<SegmentationNaive3DCloseCurve>( *segmenter3DCloseCurve, board );
	}
	else
	{
	  computer3D = new Naive3SegmentComputer();
	  segmenter3D = new SegmentationNaive3D ( digitalCurve.begin(), digitalCurve.end(), *computer3D );
	  printDSS3D<SegmentationNaive3D>( *segmenter3D, board );
	}
	
	if ( lengthFilter )
	{
	  lmst3D64CloseCurve->setLengthFilter ( lenMin/2. );
	  lmst3D64->setLengthFilter ( lenMin/2. );
	}
	if ( distanceFilter )
	{
	  lmst3D64CloseCurve->setDistanceFilter ( lenMin/2. );
	  lmst3D64->setDistanceFilter ( lenMin/2. );
	}
      }
      
      template <class TMyCurve, class TSegmentComp>
      bool testCover ( TMyCurve const & curve, TSegmentComp & dssSegments )
      {
	typename TSegmentComp::SegmentComputerIterator lastDSS = dssSegments.end();
	
	typename TMyCurve::const_iterator it = curve.begin();
	for (; it != curve.end(); ++it)
	{
	  bool flag = false;
	  typename TSegmentComp::SegmentComputerIterator DSS = dssSegments.begin();
	  for ( ; DSS != lastDSS; ++DSS )
	  {
	    if ( DSS->isInDSS ( *it ) )
	      flag = true;
	  }
	  if (!flag)
	    return false;
	}
	return true;
      }
      
      template <class TMyCurve, class TSegmentComp>
      void repairCover ( TMyCurve const & curve, TSegmentComp & dssSegments )
      {
	dssSegments.setMode("First");
	if ( !testCover ( curve, dssSegments ) )
	{
	  dssSegments.setMode("Last");
	  if ( !testCover ( curve, dssSegments ) )
	  {
	    dssSegments.setMode("MostCentered");
	    if ( !testCover ( curve, dssSegments ) )
	      throw "Curve not covered!";
	  }
	}
      }
     
      // I am not sure if we should consider all projections as is written in a draft.
      void findShortest2DMaximalSegment ()
      {
	remove2DSegements();
	lenMin = UINT_MAX;
	SimpleConstRange<typename MyProjection::const_iterator> t_RangeXY ( oXY.begin(), oXY.end() );
	SimpleConstRange<typename MyProjection::const_iterator> t_RangeXZ ( oXZ.begin(), oXZ.end() );
	SimpleConstRange<typename MyProjection::const_iterator> t_RangeYZ ( oYZ.begin(), oYZ.end() );
	
	//move to the function
	if ( isCloseCurve< MyProjection, Z2i::Point > ( oXY ) )
	{
	  segmenterXYCloseCurve = new SegmentationCloseCurve ( t_RangeXY.c(), --t_RangeXY.c(), SegmentComputerCloseCurve() );
	  repairCover<MyProjection, SegmentationCloseCurve>( oXY, *segmenterXYCloseCurve );
	  
	  for ( SegmentationCloseCurve::SegmentComputerIterator it = segmenterXYCloseCurve->begin(); it != segmenterXYCloseCurve->end(); ++it )
	    if ( std::distance ( (*it).begin(), (*it).end() ) < lenMin )
	      lenMin = std::distance ( (*it).begin(), (*it).end() );
	}
	else
	{
	  segmenterXY = new Segmentation ( oXY.begin(), oXY.end(), SegmentComputer() );
	  for ( Segmentation::SegmentComputerIterator it = segmenterXY->begin(); it != segmenterXY->end(); ++it )
	    if ( std::distance ( (*it).begin(), (*it).end() ) < lenMin )
	      lenMin = std::distance ( (*it).begin(), (*it).end() );
	}
	
	if ( isCloseCurve< MyProjection, Z2i::Point > ( oXZ ) )
	{
	  segmenterXZCloseCurve = new SegmentationCloseCurve ( t_RangeXZ.c(), --t_RangeXZ.c(), SegmentComputerCloseCurve() );
	  repairCover<MyProjection, SegmentationCloseCurve>( oXZ, *segmenterXZCloseCurve );
	  for ( SegmentationCloseCurve::SegmentComputerIterator it = segmenterXZCloseCurve->begin(); it != segmenterXZCloseCurve->end(); ++it )
	    if ( std::distance ( (*it).begin(), (*it).end() ) < lenMin )
	      lenMin = std::distance ( (*it).begin(), (*it).end() );
	}
	else
	{
	  segmenterXZ = new Segmentation ( oXZ.begin(), oXZ.end(), SegmentComputer() );
	  for ( Segmentation::SegmentComputerIterator it = segmenterXZ->begin(); it != segmenterXZ->end(); ++it )
	    if ( std::distance ( (*it).begin(), (*it).end() ) < lenMin )
	      lenMin = std::distance ( (*it).begin(), (*it).end() );
	}
	if ( isCloseCurve< MyProjection, Z2i::Point > ( oYZ ) )
	{
	  segmenterYZCloseCurve = new SegmentationCloseCurve ( t_RangeYZ.c(), --t_RangeYZ.c(), SegmentComputerCloseCurve() );
	  repairCover<MyProjection, SegmentationCloseCurve>( oYZ, *segmenterYZCloseCurve );
	  for ( SegmentationCloseCurve::SegmentComputerIterator it = segmenterYZCloseCurve->begin(); it != segmenterYZCloseCurve->end(); ++it )
	    if ( std::distance ( (*it).begin(), (*it).end() ) < lenMin )
	      lenMin = std::distance ( (*it).begin(), (*it).end() );
	}
	else
	{
	  segmenterYZ = new Segmentation ( oYZ.begin(), oYZ.end(), SegmentComputer() );
	  for ( Segmentation::SegmentComputerIterator it = segmenterYZ->begin(); it != segmenterYZ->end(); ++it )
	    if ( std::distance ( (*it).begin(), (*it).end() ) < lenMin )
	      lenMin = std::distance ( (*it).begin(), (*it).end() );
	}
      }
      
      //! \todo cleanup
      void _test ( double & mse2d, double & mse3d, double & mae2d, double & mae3d, double &mse2DVar, double &mse3DVar )
      {
	std::stringstream tr;
	tr << "HelixMST3D_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream curve3D ( tr.str() );
	tr.str("");
	tr << "HelixMST2D_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream projec ( tr.str() );
	tr.str("");
	tr << "HelixGroundTruth_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream ground ( tr.str() );
	
	
	Z3i::RealVector tanGroundTruth;
	Z3i::RealVector tanMST2D;
	
	std::vector<double> mse2dPartial;
	std::vector<double> mse3dPartial;
	
	bool flag = true;
	unsigned int axis = -1;
	unsigned int counter = 0;
	for ( double i = mint; i < maxt; i += 0.05 )
	{
	  Z3i::RealPoint realPoint = rotCurve->x ( i );
	  typename MyDigitalCurve::const_iterator it = FindClosestPoint < MyDigitalCurve > ( realPoint, digitalCurve );
	  
	  std::cout << "Point: " << *it << std::endl;
	  
	  //omit ends for open curve -- they never have a good approximation
	  if ( !( isCloseCurve<MyDigitalCurve, Z3i::Point>( digitalCurve ) ) && ( it == digitalCurve.begin() || it == digitalCurve.end() ) )
	    continue;
	  
	  axis = findMainAxis<MyRotatedCurve> ( *rotCurve, i );
	  Z2i::RealVector v0, v1;
	  if ( axis == 0 )
	  {
	    if ( isCloseCurve<MyProjection, Z2i::Point> ( oXY ) )
	    {
	      lmst64CloseCurve->attach ( *segmenterXYCloseCurve );
	      v0 = lmst64CloseCurve->eval ( myProjXY ( *it ) );
	    }
	    else
	    {
	      lmst64->attach ( *segmenterXY );
	      v0 = lmst64->eval ( myProjXY ( *it ) );
	    }
	    if ( isCloseCurve<MyProjection, Z2i::Point> ( oXZ ) )
	    {
	      lmst64CloseCurve->attach ( *segmenterXZCloseCurve );
	      v1 = lmst64CloseCurve->eval ( myProjXZ ( *it ) );
	    }
	    else
	    {
	      lmst64->attach ( *segmenterXZ );
	      v1 = lmst64->eval ( myProjXZ ( *it ) );
	    }
	    if ( v1[0] == 0 || ( v0[1] == 0 && v1[1] == 0 ) )
	    {
	      tanMST2D[0] = v0[1];
	      tanMST2D[1] = v0[0];
	      tanMST2D[2] = v1[0];
	    }
	    else
	    {
	      if ( v0[0] == 0 )
	      {
		tanMST2D[0] = v1[1];
		tanMST2D[1] = 0;
		tanMST2D[2] = v1[0];
	      }
	      else
	      {
		tanMST2D[0] = v1[1] * v0[1];
		tanMST2D[1] = v1[1] * v0[0];
		tanMST2D[2] = v0[1] * v1[0];
	      }
	    }
	    if ( tanMST2D.norm() != 0 )
	      tanMST2D = tanMST2D.getNormalized();
	  }
	  else if ( axis == 1 )
	  {
	    if ( isCloseCurve<MyProjection, Z2i::Point> ( oXY ) )
	    {
	      lmst64CloseCurve->attach ( *segmenterXYCloseCurve );
	      v0 = lmst64CloseCurve->eval ( myProjXY ( *it ) );
	    }
	    else
	    {
	      lmst64->attach ( *segmenterXY );
	      v0 = lmst64->eval ( myProjXY ( *it ) );
	    }
	    if ( isCloseCurve<MyProjection, Z2i::Point> ( oYZ ) )
	    {
	      lmst64CloseCurve->attach ( *segmenterYZCloseCurve );
	      v1 = lmst64CloseCurve->eval ( myProjYZ ( *it ) );
	    }
	    else
	    {
	      lmst64->attach ( *segmenterYZ );
	      v1 = lmst64->eval ( myProjYZ ( *it ) );
	    }
	    
	    if ( v0[1] == 0 || ( v1[1] == 0 && v0[0] == 0 ) )
	    {
	      tanMST2D[0] = v0[1];
	      tanMST2D[1] = v1[1];
	      tanMST2D[2] = v1[0];
	    }
	    else
	    {
	      if ( v1[0] == 0 )
	      {
		tanMST2D[0] = v0[1];
		tanMST2D[1] = v0[0];
		tanMST2D[2] = 0;
	      }
	      else
	      {
		tanMST2D[0] = v1[1] * v0[1];
		tanMST2D[1] = v1[1] * v0[0];
		tanMST2D[2] = v0[0] * v1[0];
	      }
	    }
	    if ( tanMST2D.norm() != 0 )
	      tanMST2D = tanMST2D.getNormalized();
	  }
	  else if ( axis == 2 )
	  {
	    if ( isCloseCurve<MyProjection, Z2i::Point> ( oXZ ) )
	    {
	      lmst64CloseCurve->attach ( *segmenterXZCloseCurve );
	      v0 = lmst64CloseCurve->eval ( myProjXZ ( *it ) );
	    }
	    else
	    {
	      lmst64->attach ( *segmenterXZ );
	      v0 = lmst64->eval ( myProjXZ ( *it ) );
	    }
	    if ( isCloseCurve<MyProjection, Z2i::Point> ( oYZ ) )
	    {
	      lmst64CloseCurve->attach ( *segmenterYZCloseCurve );
	      v1 = lmst64CloseCurve->eval ( myProjYZ ( *it ) );
	    }
	    else
	    {
	      lmst64->attach ( *segmenterYZ );
	      v1 = lmst64->eval ( myProjYZ ( *it ) );
	    }
	    
	    if ( v0[1] == 0 || ( v0[0] == 0 && v1[0] == 0 ) )
	    {
	      tanMST2D[0] = v0[1];
	      tanMST2D[1] = v1[1];
	      tanMST2D[2] = v1[0];
	    }
	    else
	    {
	      if ( v1[1] == 0 )
	      {
		tanMST2D[0] = v0[1];
		tanMST2D[1] = 0;
		tanMST2D[2]= v0[0];
	      }
	      else
	      {
		tanMST2D[0] = v0[1] * v1[0];
		tanMST2D[1] = v1[1] * v0[0];
		tanMST2D[2] = v0[0] * v1[0];
	      }
	    }
	    
	    if ( tanMST2D.norm() != 0 )
	      tanMST2D = tanMST2D.getNormalized();
	  }
	  
	  tanGroundTruth = rotCurve->xp ( i );
	  tanGroundTruth = tanGroundTruth.getNormalized();
	  counter++;
	  if ( tanMST2D[0] < 0 && tanGroundTruth[0] > 0 )
	    tanMST2D[0] = -tanMST2D[0];
	  if ( tanMST2D[1] < 0 && tanGroundTruth[1] > 0 )
	    tanMST2D[1] = -tanMST2D[1];
	  if ( tanMST2D[2] < 0 && tanGroundTruth[2] > 0 )
	    tanMST2D[2] = -tanMST2D[2];
	  
	  if ( tanMST2D[0] > 0 && tanGroundTruth[0] < 0 )
	    tanMST2D[0] = -tanMST2D[0];
	  if ( tanMST2D[1] > 0 && tanGroundTruth[1] < 0 )
	    tanMST2D[1] = -tanMST2D[1];
	  if ( tanMST2D[2] > 0 && tanGroundTruth[2] < 0 )
	    tanMST2D[2] = -tanMST2D[2];
	  
	  Z3i::RealPoint tangent3D;
	  if ( isCloseCurve<MyDigitalCurve, Z3i::Point>( digitalCurve ) )
	    tangent3D = lmst3D64CloseCurve->eval(*it);
	  else 
	    tangent3D = lmst3D64->eval(*it);
	  
	  if ( tangent3D.norm() != 0 )
	    tangent3D = tangent3D.getNormalized();
	  
	  if ( tangent3D[0] < 0 && tanGroundTruth[0] > 0 )
	    tangent3D[0] = -tangent3D[0];
	  if ( tangent3D[1] < 0 && tanGroundTruth[1] > 0 )
	    tangent3D[1] = -tangent3D[1];
	  if ( tangent3D[2] < 0 && tanGroundTruth[2] > 0 )
	    tangent3D[2] = -tangent3D[2];
	  
	  if ( tangent3D[0] > 0 && tanGroundTruth[0] < 0 )
	    tangent3D[0] = -tangent3D[0];
	  if ( tangent3D[1] > 0 && tanGroundTruth[1] < 0 )
	    tangent3D[1] = -tangent3D[1];
	  if ( tangent3D[2] > 0 && tanGroundTruth[2] < 0 )
	    tangent3D[2] = -tangent3D[2];
	  
	  double mse3DTmp = pow ( tangent3D[0] - tanGroundTruth[0], 2 ) + pow ( tangent3D[1] - tanGroundTruth[1], 2 ) + pow (tangent3D[2] - tanGroundTruth[2], 2 );
	  mse3d +=  mse3DTmp;
	  mse3dPartial.push_back ( mse3DTmp );
	  double mse2DTmp = pow ( tanMST2D[0] - tanGroundTruth[0], 2 ) + pow ( tanMST2D[1] - tanGroundTruth[1], 2 ) + pow ( tanMST2D[2] - tanGroundTruth[2], 2 );
	  mse2d +=  mse2DTmp;
	  mse2dPartial.push_back ( mse2DTmp );
	  mae3d +=  fabs ( tangent3D[0] - tanGroundTruth[0] ) + fabs ( tangent3D[1] - tanGroundTruth[1] ) + fabs (tangent3D[2] - tanGroundTruth[2] );
	  mae2d +=  fabs ( tanMST2D[0] - tanGroundTruth[0] ) + fabs ( tanMST2D[1] - tanGroundTruth[1] ) + fabs ( tanMST2D[2] - tanGroundTruth[2] );
	  curve3D << i << "," << tangent3D[0] << "," << tangent3D[1] << "," << tangent3D[2] << std::endl;
	  projec << i << "," << tanMST2D[0] << "," << tanMST2D[1] << "," << tanMST2D[2] << std::endl;
	  ground << i << "," << tanGroundTruth[0] << "," << tanGroundTruth[1] << "," << tanGroundTruth[2] << std::endl;
	}
	
	mse2d /= (double)counter;
	mse3d /= (double)counter;
	mae2d /= (double)counter;
	mae3d /= (double)counter;
	
	assert ( counter == mse2dPartial.size() );
	for ( unsigned int i = 0; i < mse2dPartial.size(); i++ )
	{
	  mse2DVar += pow ( mse2dPartial[i] - mse2d, 2 );
	  mse3DVar += pow ( mse3dPartial[i] - mse2d, 2 );
	}
	mse2DVar /= mse2dPartial.size();
	mse3DVar /= mse3dPartial.size();
	
	curve3D.close();
	projec.close();
	ground.close();
      }
      
public:
  TestRunner ( Board3D<> & p_board) : board(p_board)
  {
    radius1 = radius2 = depth = angle = mse2d = mae3d = mae2d = mae3d = 0.;
    lengthFilter = false;
    distanceFilter = false;
    realCurve  = nullptr;
    trans = nullptr;
    inverse = nullptr;
    rotCurve = nullptr;
    segmenterXYCloseCurve = nullptr;
    segmenterXY = nullptr;
    segmenterXZCloseCurve = nullptr;
    segmenterXZ = nullptr;
    segmenterYZCloseCurve = nullptr;
    segmenterYZ = nullptr;
    computer3DCloseCurve = nullptr;
    segmenter3DCloseCurve = nullptr;
    computer3D = nullptr;
    segmenter3D = nullptr;
    lmst3D64CloseCurve = nullptr;
    lmst3D64 = nullptr;
    lmst64CloseCurve = nullptr;
    lmst64 = nullptr;
    std::vector<DGtal::Dimension> v1,v2,v3;
    v1.push_back(0);
    v1.push_back(1);
    v2.push_back(0);
    v2.push_back(2);
    v3.push_back(1);
    v3.push_back(2);
    myProjXY.init(v1.begin(),v1.end());
    myProjXZ.init(v2.begin(),v2.end());
    myProjYZ.init(v3.begin(),v3.end());
    lenMin = UINT_MAX;
  }
  
  void setCurveParams ( float p_radius1, float p_radius2, float p_depth )
  {
    radius1 = p_radius1;
    radius2 = p_radius2;
    depth = p_depth;
    
    delete realCurve;
    realCurve = new TCurve ( radius1, radius2, depth );
  }
  
  void setTransformationParams ( double p_angle, Z3i::RealVector p_axis )
  {
    angle = p_angle;
    axis = p_axis;
    
    delete trans;
    delete inverse;
    delete rotCurve;
    
    trans = new ForwardTrans ( RealPoint ( 0, 0, 0 ), axis, angle, RealVector (0,0,0) );
    inverse = new BackwardTrans ( RealPoint ( 0, 0, 0 ), axis, angle, RealVector (0,0,0) );
    rotCurve = new MyRotatedCurve ( *realCurve, *trans, *inverse );
  }
  
  void setLengthFilter (bool filter)
  {
    lengthFilter = filter;
  }
  
  void setDistanceFilter (bool filter)
  {
    distanceFilter = filter;
  }
  
  void digitizeCurve ( double tmin, double tmax, bool changeAxis, bool & status )
  {
    isInitialized();
    mint = tmin;
    maxt = tmax;
    digitalCurve.clear();
    Digitizer digitizer;
    digitizer.attach ( *rotCurve );
    digitizer.init ( digitalCurve );
    try
    {
      digitizer.digitize ( mint, maxt, 0.00001 );
    } catch (...)
    {
      status = false;
      if ( maxt < ( tmax / 2. ) || digitalCurve.size() < 5 )
	throw "Curve too short.";
    }
    generateProjections();
    exportPINK<MyDigitalCurve> ( radius1, radius2, depth, digitalCurve );
    print3D<MyDigitalCurve> ( digitalCurve, board );
  }
  
  void runTest ( double & mse2d, double & mse3d, double & mae2d, double & mae3d, double &mse2dVar, double &mse3DVar )
  {
    generateProjections();
    findShortest2DMaximalSegment();
    initEstimators();
    _test ( mse2d, mse3d, mae2d, mae3d, mse2dVar, mse3DVar );
  }
  
  ~TestRunner()
  {
    delete trans;
    delete inverse;
    delete rotCurve;
    remove2DSegements();
    remove3DSegements();
    removeEstimators();
  }
};

int main ( int argc, char ** argv )
{
  if ( argc < 5 )
    throw std::runtime_error ( "Too few arguments!" );
//      feenableexcept ( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW );
  std::ofstream mse ( "mse_ell_helix.csv" );
  mse << "ln(h)" << "," << "ln(mse2d)" << "," << "ln(mse3d)" << "," << "h" << "," << "mse2d" << "," << "mse3d" << "," << "mae2d" << "," << "mae3d" << "," << "mse2dVar" << "," << "mse3dVar" << std::endl;
  double mse2d, mse3d, mae2d, mae3d, mse2dVar, mse3DVar;
  Board3D<> board;
  TestRunner< EllipticHelix < Z3i::Space > > testRunner ( board );
  for ( double i = 1.;  i < std::atof ( argv[5] ); i += std::atof ( argv[4] ) )
  {
    mse2d = mse3d = mae2d = mae3d = mse2dVar = mse3DVar = 0.;
    bool status = true;
    board.clear();
    std::stringstream tr;
    tr << std::setprecision(3) << std::atof ( argv[1] ) * i << "_" << std::setprecision(3) << std::atof ( argv[2] ) * i << "_" << std::setprecision(3) << std::atof ( argv[3] ) * i;
//     testRunner.setLengthFilter ( true );
//     testRunner.setDistanceFilter ( true );
    testRunner.setCurveParams ( std::atof ( argv[1] ) * i, std::atof ( argv[2] ) * i,std::atof ( argv[3] ) * i );
    testRunner.setTransformationParams ( std::atof(argv[6]), Z3i::RealVector ( 1, 0, 1 )  );
    testRunner.digitizeCurve ( 0.1, 2. * M_PI - 0.1, true, status );
    testRunner.runTest ( mse2d, mse3d, mae2d, mae3d, mse2dVar, mse3DVar );
    board.saveOBJ( tr.str().c_str() );
    mse << std::log ( i ) << "," << std::log (mse2d) << "," << std::log ( mse3d ) << "," << i << "," << mse2d << "," << mse3d << "," << mae2d << "," << mae3d
    << "," << mse2dVar << "," << mse3DVar << "," << status << std::endl;;
  }
  std::cout << "Done." << std::endl;
  mse.close();
  return 0;
}