#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <vector>
#include<DGtal/helpers/StdDefs.h>
#include "DGtal/curves/ParametricCurveDigitizer3D.h"
#include "DGtal/curves/parametric/EllipticHelix.h"
#include <DGtal/geometry/curves/SaturatedSegmentation.h>
#include "DGtal/curves/parametric/DecoratorCurveTransformation.h"
#include "DGtal/images/RigidTransformation3D.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/Naive3DDSSComputer.h"
#include "DGtal/geometry/curves/LambdaMST2D.h"
#include "DGtal/geometry/curves/LambdaMST3D.h"
#include "DGtal/geometry/curves/LambdaMST3DBy2D.h"
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
    board.setFillColor ( Color ( 0, 0, 255, 100 ) );
    board << SetMode3D(digitalCurve[i].className(), "");
    board << digitalCurve[i];
    //     board << SetMode3D(digitalCurve[i].className(), "PavingWired");
    //     board << digitalCurve[i];
  }
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

template <typename TPoint>
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

template <typename TSegmentComp>
void printDSS3D ( TSegmentComp & segmenter, Board3D<> & board )
{
  for ( typename TSegmentComp::SegmentComputerIterator it = segmenter.begin() ; it != segmenter.end(); ++it )
  {
    board << SetMode3D((*it).className(), "BoundingBox");
    board.setFillColor ( Color ( 255, 255, 255, 255 ) );
    board << *it;
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
      Naive3SegmentComputerCloseCurve * computer3DCloseCurve;
      SegmentationNaive3DCloseCurve * segmenter3DCloseCurve;
      Naive3SegmentComputer * computer3D;
      SegmentationNaive3D * segmenter3D;
      LambdaMST3D < SegmentationNaive3DCloseCurve > * lmst3D64CloseCurve;
      LambdaMST3D < SegmentationNaive3D > * lmst3D64;
      //!\todo add version with circulator
      LambdaMST3DBy2D < typename MyDigitalCurve::const_iterator > * lmst64By2D;
      LambdaMST3DBy2D < typename SimpleConstRange < typename MyDigitalCurve::const_iterator >::ConstCirculator > * lmst64By2DCloseCurve;
      unsigned int lenMin;
      Board3D<> & board;
      Board2D board2d;
      
      double mint, maxt;
      
      void isInitialized()
      {
	if ( realCurve == nullptr || trans == nullptr || inverse == nullptr || rotCurve == nullptr  )
	  throw "Not initialized! First call; setCurveParams; setTransformationParams.";
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
	delete lmst64By2D;
	lmst64By2D = nullptr;
	delete lmst64By2DCloseCurve;
	lmst64By2DCloseCurve = nullptr;
	delete lmst3D64CloseCurve,
	lmst3D64CloseCurve = nullptr,
	delete lmst3D64;
	lmst3D64 = nullptr;
      }
      
      void initEstimators()
      {
	remove3DSegements();
	removeEstimators();
	
	lmst64By2D = new LambdaMST3DBy2D < typename MyDigitalCurve::const_iterator >();
	lmst64By2DCloseCurve = new LambdaMST3DBy2D < typename SimpleConstRange < typename MyDigitalCurve::const_iterator >::ConstCirculator >();
	
	computer3D = new Naive3SegmentComputer();
	segmenter3D = new SegmentationNaive3D ( digitalCurve.begin(), digitalCurve.end(), *computer3D );
	lmst3D64CloseCurve = new LambdaMST3D < SegmentationNaive3DCloseCurve >();
	lmst3D64 = new LambdaMST3D < SegmentationNaive3D >();
	lmst3D64->init ( digitalCurve.begin(), digitalCurve.end() );
	lmst3D64->attach ( *segmenter3D );
	if ( isCloseCurve<MyDigitalCurve, Z3i::Point> ( digitalCurve ) )
	{
	  SimpleConstRange<typename MyDigitalCurve::const_iterator> t_Range(digitalCurve.begin(), digitalCurve.end());
	  computer3DCloseCurve = new Naive3SegmentComputerCloseCurve();
	  segmenter3DCloseCurve = new SegmentationNaive3DCloseCurve ( t_Range.c(), --t_Range.c(), *computer3DCloseCurve );
	  segmenter3DCloseCurve->setMode ( "MostCentered" );
	  repairCover<MyDigitalCurve, SegmentationNaive3DCloseCurve>( digitalCurve, *segmenter3DCloseCurve );
	  lmst3D64CloseCurve->init ( t_Range.c(), --t_Range.c() );
	  lmst3D64CloseCurve->attach ( *segmenter3DCloseCurve );
	  // 	  printDSS3D<SegmentationNaive3DCloseCurve>( *segmenter3DCloseCurve, board );
	}
	else
	{
	  computer3D = new Naive3SegmentComputer();
	  segmenter3D = new SegmentationNaive3D ( digitalCurve.begin(), digitalCurve.end(), *computer3D );
	  // 	  printDSS3D<SegmentationNaive3D>( *segmenter3D, board );
	}
      }
      
      void findShortestValid2DMinimalSegment()
      {
	unsigned int lenXY = 0; unsigned int lenXZ = 0; unsigned int lenYZ = 0;
	lenMin = UINT_MAX;
	if ( isCloseCurve<MyDigitalCurve, Z3i::Point> ( digitalCurve ) )
	{
	  typename SegmentationNaive3DCloseCurve::SegmentComputerIterator it = segmenter3DCloseCurve->begin();
	  for (; it != segmenter3DCloseCurve->end(); ++it )
	  {
	    (*it).get2DSegmentsLength ( lenXY, lenXZ, lenYZ );
	    if ( ( lenXY > lenYZ && lenXZ > lenYZ ) && lenMin > std::min ( lenXY, lenXZ ) )
	      lenMin = std::min ( lenXY, lenXZ );
	    else if ( ( lenYZ > lenXZ && lenXY > lenXZ ) && lenMin > std::min ( lenXY, lenYZ ) )
	      lenMin = std::min ( lenXY, lenYZ );
	    else if ( ( lenYZ > lenXY && lenXZ > lenXY ) && lenMin > std::min ( lenXZ, lenYZ ) )
	      lenMin = std::min ( lenXZ, lenYZ );
	  }
	}
	else
	{
	  typename SegmentationNaive3D::SegmentComputerIterator it = segmenter3D->begin();
	  for (; it != segmenter3D->end(); ++it )
	  {
	    (*it).get2DSegmentsLength ( lenXY, lenXZ, lenYZ );
	    if ( ( lenXY > lenYZ && lenXZ > lenYZ ) && lenMin > std::min ( lenXY, lenXZ ) )
	      lenMin = std::min ( lenXY, lenXZ );
	    else if ( ( lenYZ > lenXZ && lenXY > lenXZ ) && lenMin > std::min ( lenXY, lenYZ ) )
	      lenMin = std::min ( lenXY, lenYZ );
	    else if ( ( lenYZ > lenXY && lenXZ > lenXY ) && lenMin > std::min ( lenXZ, lenYZ ) )
	      lenMin = std::min ( lenXZ, lenYZ );
	  }
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
	if ( !testCover ( curve, dssSegments ) )
	{
	  dssSegments.setMode("First");
	  if ( !testCover ( curve, dssSegments ) )
	  {
	    dssSegments.setMode("Last");
	    if ( !testCover ( curve, dssSegments ) )
	      throw "Curve not covered!";
	  }
	}
      }
      
      void initEstimatorsBy2D ()
      {
	SimpleConstRange < typename MyDigitalCurve::const_iterator > t_Range ( digitalCurve.begin(), digitalCurve.end() );
	lmst64By2DCloseCurve->init ( t_Range.c(), --t_Range.c() );
	lmst64By2D->init ( digitalCurve.begin(), digitalCurve.end() );
      }
      
      void _test ( double & mse2d, double & mse3d, double & mse3dF,
		   double & mae2d, double & mae3d, double & mae3dF,
		   double & mse2DVar, double & mse3DVar, double & mse3DVarF )
      {
	std::stringstream tr;
	tr << "HelixMST3D_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream curve3D ( tr.str() );
	tr.str("");
	tr << "HelixMST3DFiltered_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream curve3DFiltered ( tr.str() );
	tr.str("");
	tr << "HelixMST2D_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream projec ( tr.str() );
	tr.str("");
	tr << "HelixGroundTruth_" << std::showpoint << std::setprecision(3) << radius1 << "_" << std::showpoint << std::setprecision(3) << radius2 << "_" << std::showpoint << std::setprecision(3) << depth << ".csv";
	std::ofstream ground ( tr.str() );
	
	Z3i::RealVector tanGroundTruth;
	Z3i::RealPoint tangent3D;
	Z3i::RealPoint tangent3DF;
	PointVector<3, double> tanMST2D;
	
	std::vector<double> mse2dPartial;
	std::vector<double> mse3dPartial;
	std::vector<double> mse3dFPartial;
	std::vector<Z3i::Point> discardedSet;
	
	bool flag = true;
	unsigned char axis = -1, axisPrev = -1;
	unsigned int counter = 0;
	for ( double i = mint; i < maxt; i += 0.05 )
	{
	  Z3i::RealPoint realPoint = rotCurve->x ( i );
	  tanGroundTruth = rotCurve->xp ( i );
	  tanGroundTruth = tanGroundTruth.getNormalized();
	  typename MyDigitalCurve::const_iterator it = FindClosestPoint < MyDigitalCurve > ( realPoint, digitalCurve );
	  
	  std::cout << "Point: " << *it << std::endl;
	  
	  //omit ends for open curve -- they never have a good approximation
	  if ( !( isCloseCurve<MyDigitalCurve, Z3i::Point>( digitalCurve ) ) && ( it == digitalCurve.begin() || it == digitalCurve.end() ) )
	    continue;
	  
	  //Do not consider points which are too far from real curve.
	  if ( distance ( RealPoint ( *it ), realPoint  ) >= std::sqrt ( 3 ) / 2. )
	    continue;
	  
	  //check if director changed and manage discarded set.
	  axisPrev = axis;
	  axis = findMainAxis<MyRotatedCurve>(*rotCurve, i);
	  if (  axis != axisPrev )
	  {
	    if ( std::find ( discardedSet.begin(), discardedSet.end(), *it ) == discardedSet.end() )
	      discardedSet.push_back ( *it );
	  }
	  
	  if ( isCloseCurve<MyDigitalCurve, Z3i::Point>( digitalCurve ) )
	  {
	    // check if points are included in DSSs with points in discarded set.
	    if ( InDisartedSet < std::vector<Z3i::Point>::const_iterator, SegmentationNaive3DCloseCurve, Z3i::Point  > ( discardedSet.cbegin(), discardedSet.cend(), segmenter3DCloseCurve, *it ) )
	      continue;
	    
	    lmst3D64CloseCurve->setLengthFilter ( 0 );
	    lmst3D64CloseCurve->setDistanceFilter ( -1 );
	    tangent3D = lmst3D64CloseCurve->eval(*it);
	    lmst3D64CloseCurve->setLengthFilter ( lenMin/2. );
	    lmst3D64CloseCurve->setDistanceFilter ( lenMin/2. );
	    tangent3DF = lmst3D64CloseCurve->eval(*it);
	    tanMST2D = lmst64By2DCloseCurve->eval(*it);
	  }
	  else
	  {
	    // check if points are included in DSSs with points in discarded set.
	    if ( InDisartedSet < std::vector<Z3i::Point>::const_iterator, SegmentationNaive3D, Z3i::Point  > ( discardedSet.cbegin(), discardedSet.cend(), segmenter3D, *it ) )
	      continue;
	    
	    lmst3D64->setLengthFilter ( 0 );
	    lmst3D64->setDistanceFilter ( -1 );
	    tangent3D = lmst3D64->eval(*it);
	    lmst3D64->setLengthFilter ( lenMin/2. );
	    lmst3D64->setDistanceFilter ( lenMin/2. );
	    tangent3DF = lmst3D64->eval(*it);
	    tanMST2D = lmst64By2D->eval(*it);
	  }
	  if ( tangent3D.norm() != 0 )
	    tangent3D = tangent3D.getNormalized();
	  
	  if ( tanMST2D.norm() != 0 )
	    tanMST2D = tanMST2D.getNormalized();
	  
	  counter++;
	  double sinSquareReal = pow ( realPoint.crossProduct ( tanGroundTruth ).norm() / ( realPoint.norm() * tanGroundTruth.norm() ), 2 );
	  double sinSquare3DMST = pow ( realPoint.crossProduct ( tangent3D ).norm() / ( realPoint.norm() * tangent3D.norm() ), 2 );
	  double sinSquare3DMSTF = pow ( realPoint.crossProduct ( tangent3DF ).norm() / ( realPoint.norm() * tangent3DF.norm() ), 2 );
	  double sinSquare3DMSTBy2D = pow ( realPoint.crossProduct ( tanMST2D ).norm() / ( realPoint.norm() * tanMST2D.norm() ), 2 );
	  
	  double mse3DTmp = pow ( sinSquare3DMST - sinSquareReal, 2 );
	  mse3d +=  mse3DTmp;
	  mse3dPartial.push_back ( mse3DTmp );
	  double mse3DTmpF = pow ( sinSquare3DMSTF - sinSquareReal, 2 );
	  mse3dF +=  mse3DTmpF;
	  mse3dFPartial.push_back ( mse3DTmpF );
	  double mse2DTmp = pow ( sinSquare3DMSTBy2D - sinSquareReal, 2 );
	  mse2d +=  mse2DTmp;
	  mse2dPartial.push_back ( mse2DTmp );
	  mae3d +=  fabs ( sinSquare3DMST - sinSquareReal );
	  mae3dF +=  fabs ( sinSquare3DMSTF - sinSquareReal );
	  mae2d +=  fabs ( sinSquare3DMSTBy2D - sinSquareReal );
	  
	  if ( tanMST2D.cosineSimilarity ( tanGroundTruth ) > M_PI_2 )
	    tanMST2D = -tanMST2D;
	  if ( tangent3D.cosineSimilarity ( tanGroundTruth ) > M_PI_2 )
	  {
	    tangent3D = -tangent3D;
	    tangent3DF = -tangent3DF;
	  }
	  curve3D << i << "," << tangent3D[0] << "," << tangent3D[1] << "," << tangent3D[2] << std::endl;
	  curve3DFiltered << i << "," << tangent3DF[0] << "," << tangent3DF[1] << "," << tangent3DF[2] << std::endl;
	  projec << i << "," << tanMST2D[0] << "," << tanMST2D[1] << "," << tanMST2D[2] << std::endl;
	  ground << i << "," << tanGroundTruth[0] << "," << tanGroundTruth[1] << "," << tanGroundTruth[2] << std::endl;
	}
	
	mse2d /= (double)counter;
	mse3d /= (double)counter;
	mse3dF /= (double)counter;
	mae2d /= (double)counter;
	mae3d /= (double)counter;
	mae3dF /= (double)counter;
	
	assert ( counter == mse2dPartial.size() );
	for ( unsigned int i = 0; i < mse2dPartial.size(); i++ )
	{
	  mse2DVar += pow ( mse2dPartial[i] - mse2d, 2 );
	  mse3DVar += pow ( mse3dPartial[i] - mse3d, 2 );
	  mse3DVarF += pow ( mse3dFPartial[i] - mse3dF, 2 );
	}
	mse2DVar /= mse2dPartial.size();
	mse3DVar /= mse3dPartial.size();
	mse3DVarF /= mse3dFPartial.size();
	
	curve3D.close();
	curve3DFiltered.close();
	projec.close();
	ground.close();
      }
      
      template <class IteratorDisartedSet, class TSegmentation, class TPoint>
      bool InDisartedSet ( const IteratorDisartedSet & dbegin, const IteratorDisartedSet & dend,
			   const TSegmentation * segmenter, const TPoint & point )
      {
	typename TSegmentation::SegmentComputerIterator itt = segmenter->begin();
	bool flag = false;
	for (; itt != segmenter->end(); ++itt )
	{
	  for ( IteratorDisartedSet it = dbegin; it != dend; ++it )
	    if ( itt->isInDSS ( point ) &&  itt->isInDSS ( *it ) )
	      return true;
	}
	return false;
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
    computer3DCloseCurve = nullptr;
    segmenter3DCloseCurve = nullptr;
    computer3D = nullptr;
    segmenter3D = nullptr;
    lmst3D64CloseCurve = nullptr;
    lmst3D64 = nullptr;
    lmst64By2D = nullptr;
    std::vector<DGtal::Dimension> v1,v2,v3;
    v1.push_back(0);
    v1.push_back(1);
    v2.push_back(0);
    v2.push_back(2);
    v3.push_back(1);
    v3.push_back(2);
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
    exportPINK<MyDigitalCurve> ( radius1, radius2, depth, digitalCurve );
    //     print3D<MyDigitalCurve> ( digitalCurve, board );
  }
  
  /**
   * @param mse2d mean square error of sin^2 between estimated tangent ( LambdaMST3DBy2D ) and point of real curve.
   * @param mse3d mean square error of sin^2 between estimated tangent ( LambdaMST3D ) and point of real curve.
   * @param mse3dF mean square error of sin^2 between estimated tangent ( LambdaMST3D filtered ) and point of real curve.
   * 
   * @param mae2d mean absolute error of sin^2 between estimated tangent ( LambdaMST3DBy2D ) and point of real curve.
   * @param mae3d mean absolute error of sin^2 between estimated tangent ( LambdaMST3D ) and point of real curve.
   * @param mae3dF mean absolute error of sin^2 between estimated tangent ( LambdaMST3D filtered ) and point of real curve.
   *
   * @param mae2DVar variance of mean absolute error of sin^2 between estimated tangent ( LambdaMST3DBy2D ) and point of real curve.
   * @param mae3DVar variance of  mean absolute error of sin^2 between estimated tangent ( LambdaMST3D ) and point of real curve.
   * @param mae3DVarF variance of mean absolute error of sin^2 between estimated tangent ( LambdaMST3D filtered ) and point of real curve.
   */ 
  void runTest  ( double & mse2d, double & mse3d, double & mse3dF,
		  double & mae2d, double & mae3d, double & mae3dF,
		  double & mse2DVar, double & mse3DVar, double & mse3DVarF )
  {
    initEstimators();
    initEstimatorsBy2D ();
    findShortestValid2DMinimalSegment();
    _test ( mse2d, mse3d, mse3dF, mae2d, mae3d, mae3dF, mse2DVar, mse3DVar, mse3DVarF );
  }
  
  ~TestRunner()
  {
    delete trans;
    delete inverse;
    delete rotCurve;
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
  mse << "ln(h)" << "," << "ln(mse2d)" << "," << "ln(mse3d)" << "," << "ln(mse3dFiltered)" << "," 
  << "h" << "," << "mse2d" << "," << "mse3d" << "," << "mse3dFiltered" << "," 
  << "mae2d" << "," << "mae3d" << "," << "mae3dFiltered" << ","
  << "mse2dVar" << "," << "mse3dVar" << "," << "mse3dVarFiltered" << std::endl;
  double mse2d, mse3d, mse3dF, mae2d, mae3d, mae3dF, mse2dVar, mse3DVar, mse3DVarF;
  Board3D<> board;
  TestRunner< EllipticHelix < Z3i::Space > > testRunner ( board );
  for ( double i = 1.;  i < std::atof ( argv[5] ); i += std::atof ( argv[4] ) )
  {
    mse2d = mse3d = mse3dF = mae2d = mae3d = mae3dF = mse2dVar = mse3DVar = mse3DVarF = 0.;
    bool status = true;
    board.clear();
    std::stringstream tr;
    tr << std::setprecision(3) << std::atof ( argv[1] ) * i << "_" << std::setprecision(3) << std::atof ( argv[2] ) * i << "_" << std::setprecision(3) << std::atof ( argv[3] ) * i;
    testRunner.setCurveParams ( std::atof ( argv[1] ) * i, std::atof ( argv[2] ) * i,std::atof ( argv[3] ) * i );
    testRunner.setTransformationParams ( std::atof(argv[6]), Z3i::RealVector ( 0, 1, 0 )  );
    testRunner.digitizeCurve ( 0., 2. * M_PI - 0., true, status );
    testRunner.runTest ( mse2d, mse3d, mse3dF, mae2d, mae3d, mae3dF, mse2dVar, mse3DVar, mse3DVarF );
    board.saveOBJ( tr.str().c_str() );
    mse << std::log ( i ) << "," << std::log (mse2d) << "," << std::log ( mse3d )  << "," << std::log ( mse3dF )
    << "," << i << "," << mse2d << "," << mse3d  << "," << mse3dF 
    << "," << mae2d << "," << mae3d  << "," << mae3dF
    << "," << mse2dVar << "," << mse3DVar << "," << mse3DVarF << "," << status << std::endl;
  }
  std::cout << "Done." << std::endl;
  mse.close();
  return 0;
}