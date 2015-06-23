/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file LambdaMST3DBy2D.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2015/06/16
 *
 * This file is part of the DGtal library.
 */

#if defined(LAMBDAMST3DBy2D_RECURSES)
#error Recursive header files inclusion detected in LambdaMST3DBy2D.h
#else // defined(LAMBDAMST3DBy2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LAMBDAMST3DBy2D_RECURSES

#if !defined LAMBDAMST3DBy2D_h
/** Prevents repeated inclusion of headers. */
#define LAMBDAMST3DBy2D_h

#include <algorithm>
#include <cmath>
#include <list>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/curves/LambdaMST2D.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"

namespace DGtal {

  template < typename Iterator3D, typename Functor, typename LambdaFunctor, int CONNECTIVITY = 8 >
class LambdaMST3DBy2DEstimator
{
public:
    // ----------------------- Types ------------------------------
public:
    typedef PointVector<3, double> RealVector3D;
    typedef PointVector<3, int> Point3D;
    typedef PointVector<2, int> Point2D;
    typedef PointVector<2, double> RealVector2D;
    typedef std::list<Point2D> TCurve2D;
    typedef ArithmeticalDSSComputer < typename TCurve2D::const_iterator, int, CONNECTIVITY > SegmentComputer2D;
    typedef SaturatedSegmentation<SegmentComputer2D> Segmentation2D;
    
    // ----------------------- Private types ------------------------------
private:
    typedef LambdaMST2D < Segmentation2D, LambdaFunctor > TEstimator;
    typedef typename Functor::MAIN_AXIS MAIN_AXIS;
    typedef functors::Projector < SpaceND < 2, int > > Projector2d;

  // ----------------------- Standard services ------------------------------
public:
  LambdaMST3DBy2DEstimator(): myBegin(), myEnd()
  {
    //projections
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
  }
    
    /**
     * Initialisation.
     * @param itb, begin iterator
     * @param ite, end iterator
     */
    void init ( const Iterator3D & itB, const Iterator3D & itE )
    {
      myBegin = itB;
      myEnd = itE;
    }

    /**
     * @param point to calculate A and B for it
     * @return A and B
     */
    RealVector3D eval ( const Point3D & point )
    {
      assert ( isValid() );
      Iterator3D it = std::find ( myBegin, myEnd, point );
      TCurve2D tXY, tXZ, tYZ;
    
      ExtendBack ( tXY, tXZ, tYZ, it );
      ExtendFront ( tXY, tXZ, tYZ, it );
      
      MAIN_AXIS axis = detectMainAxis ( tXY, tXZ, tYZ, point );
      if ( axis == MAIN_AXIS::X )
	return myFunctor ( MAIN_AXIS::X, Estimate2DTangent ( tXY, myProjXY ( *it ) ), Estimate2DTangent ( tXZ, myProjXZ ( *it ) ) );
      else if ( axis == MAIN_AXIS::Y )
	return myFunctor ( MAIN_AXIS::Y, Estimate2DTangent ( tXY, myProjXY ( *it ) ), Estimate2DTangent ( tYZ, myProjYZ ( *it ) ) );
      else
	return myFunctor ( MAIN_AXIS::Z, Estimate2DTangent ( tXZ, myProjXZ ( *it ) ), Estimate2DTangent ( tYZ, myProjYZ ( *it ) ) );
    }

    /**
     * @param result output iterator on the estimated quantity
     * @todo find better algorithm
     *
     * @return the estimated quantity
     * from itb till ite (excluded)
     */
    template <typename ResultType>
    void eval ( std::vector<ResultType> & result )
    {
      assert ( isValid() );
      for( Iterator3D it = myBegin; it != myEnd; ++it )
      {
	RealVector3D tan = eval ( *it );
	if ( result.size() > 0 && tan.cosineSimilarity ( result.back() ) > M_PI_2 )
	  tan = -tan;
	result.push_back ( tan );
      }
    }
    
    // ----------------------- Standard services ------------------------------
public:

   /**
   * Checks the validity/consistency of the object.
   * @return 'true' if the object is valid, 'false' otherwise.
   */
    bool isValid() const
    {
        return ( myBegin != myEnd );
    }
    
    // ------------------------- Internals ------------------------------------
protected:
  
  inline
  MAIN_AXIS detectMainAxis ( const TCurve2D & tXY, const TCurve2D & tXZ, const TCurve2D & tYZ,
			     const Point3D & point )
  {
    unsigned int rankXY = CurveRank ( tXY, myProjXY ( point ) );
    unsigned int rankXZ = CurveRank ( tXZ, myProjXZ ( point ) );
    unsigned int rankYZ = CurveRank ( tYZ, myProjYZ ( point ) );
    if ( rankXY > rankYZ && rankXZ > rankYZ )
      return MAIN_AXIS::X;
    else if ( rankXY > rankXZ && rankYZ > rankXZ )
      return MAIN_AXIS::Y;
    else
      return MAIN_AXIS::Z;
  }
  
  unsigned int CurveRank ( const TCurve2D & curve, const Point2D & point )
  {
    unsigned int rank = 0;
    Segmentation2D segmenter ( curve.begin(), curve.end(), SegmentComputer2D() );
    for ( typename Segmentation2D::SegmentComputerIterator it = segmenter.begin(); it != segmenter.end(); ++it )
      if ( it->isInDSS ( point ) )
	rank += std::distance ( it->begin(), it->end() );
      return rank;
  }
  
  void ExtendFront ( TCurve2D & curveXY, TCurve2D & curveXZ, TCurve2D & curveYZ, const Iterator3D & it )
  {
    Iterator3D front = it;
    bool status = true, xy = true, xz = true, yz = true;
    while ( status && front >= myBegin && front < myEnd )
    {
      unsigned int test = 0;
      if ( xy && myProjXY (*front) != curveXY.front() )
      {
	curveXY.push_front ( myProjXY ( *front ) );
	test++;
      }
      else
	xy = false;
      if ( xz && myProjXZ (*front) != curveXZ.front() )
      {
	curveXZ.push_front ( myProjXZ ( *front ) );
	test++;
      }
      else
	xz = false;
      if ( yz && myProjYZ (*front) != curveYZ.front() )
      {
	curveYZ.push_front ( myProjYZ ( *front ) );
	test++;
      }
      else
	yz = false;
      if ( test >= 2 )
      {
	--front;
	status = true;
      }
      else
	status = false;
    }
  }
  
  void ExtendBack ( TCurve2D & curveXY, TCurve2D & curveXZ, TCurve2D & curveYZ, const Iterator3D & it )
  {
    Iterator3D back = it; ++back;
    bool status = true, xy = true, xz = true, yz = true;
    while ( status && back > myBegin && back < myEnd )
    {
      unsigned int test = 0;
      if ( xy && myProjXY (*back) != curveXY.back() )
      {
	curveXY.push_back ( myProjXY ( *back ) );
	test++;
      }
      else
	xy = false;
      if ( xz && myProjXZ (*back) != curveXZ.back() )
      {
	curveXZ.push_back ( myProjXZ ( *back ) );
	test++;
      }
      else
	xz = false;
      if ( yz && myProjYZ (*back) != curveYZ.back() )
      {
	curveYZ.push_back ( myProjYZ ( *back ) );
	test++;
      }
      else
	yz = false;
      if ( test >= 2 )
      {
	++back;
	status = true;
      }
      else
	status = false;
    }
  }
  
  RealVector2D Estimate2DTangent ( const TCurve2D & curve, const Point2D & point )
  {
    Segmentation2D segmenter ( curve.begin(), curve.end(), SegmentComputer2D() );
    TEstimator lmst;
    lmst.attach ( segmenter );
    return lmst.eval ( point ); 
  }

    // ------------------------- Private Datas --------------------------------
private:
    Iterator3D myBegin;
    Iterator3D myEnd;
    Functor myFunctor;
    /// projectors
    Projector2d myProjXY, myProjXZ, myProjYZ;
}; // end of class LambdaMST3DBy2DEstimator


/**
* Description of class 'LambdaTangentFromDSS' <p> Aim:
* \todo -- description
*/
   class TangentFromDSS3DBy2DFunctor
   {
   public:
     // ----------------------- Types ------------------------------
     typedef PointVector<3, double> Vector3D;
     typedef PointVector<2, double> Vector2D;
     enum MAIN_AXIS {X = 0, Y = 1, Z = 2};
     
     // ----------------------- Interface --------------------------------------
     Vector3D operator() ( MAIN_AXIS mainAxis, const Vector2D & v0, const Vector2D & v1 ) const
     {
       Vector3D tangent;
       if ( mainAxis == X )
       {
	 if ( v1[0] == 0 || ( v0[1] == 0 && v1[1] == 0 ) )
	 {
	   tangent[0] = v0[1];
	   tangent[1] = v0[0];
	   tangent[2] = v1[0];
	 }
	 else
	 {
	   if ( v0[0] == 0 )
	   {
	     tangent[0] = v1[1];
	     tangent[1] = 0;
	     tangent[2] = v1[0];
	   }
	   else
	   {
	     tangent[0] = v1[1] * v0[1];
	     tangent[1] = v1[1] * v0[0];
	     tangent[2] = v0[1] * v1[0];
	   }
	 }
       }
       else if ( mainAxis == Y )
       {
	 if ( v0[1] == 0 || ( v1[1] == 0 && v0[0] == 0 ) )
	 {
	   tangent[0] = v0[1];
	   tangent[1] = v1[1];
	   tangent[2] = v1[0];
	 }
	 else
	 {
	   if ( v1[0] == 0 )
	   {
	     tangent[0] = v0[1];
	     tangent[1] = v0[0];
	     tangent[2] = 0;
	   }
	   else
	   {
	     tangent[0] = v1[1] * v0[1];
	     tangent[1] = v1[1] * v0[0];
	     tangent[2] = v0[0] * v1[0];
	   }
	 }
       }
       else
       {
	 if ( v0[1] == 0 || ( v0[0] == 0 && v1[0] == 0 ) )
	 {
	   tangent[0] = v0[1];
	   tangent[1] = v1[1];
	   tangent[2] = v1[0];
	 }
	 else
	 {
	   if ( v1[1] == 0 )
	   {
	     tangent[0] = v0[1];
	     tangent[1] = 0;
	     tangent[2]= v0[0];
	   }
	   else
	   {
	     tangent[0] = v0[1] * v1[0];
	     tangent[1] = v1[1] * v0[0];
	     tangent[2] = v0[0] * v1[0];
	   }
	 }
       }
       return tangent;
     }
   };

//-------------------------------------------------------------------------------------------

// Template class LambdaMST3D
/**
 * \brief Aim: Simplify creation of Lambda MST tangent estimator.
 *
 */
template < 
typename Iterator3D, typename LambdaFunctor = functors::Lambda64Function, int CONNECTIVITY = 8 >
class LambdaMST3DBy2D:
public LambdaMST3DBy2DEstimator < Iterator3D, TangentFromDSS3DBy2DFunctor, LambdaFunctor, CONNECTIVITY >
  {
    typedef LambdaMST3DBy2DEstimator < Iterator3D, TangentFromDSS3DBy2DFunctor, LambdaFunctor, CONNECTIVITY > Super;
    
  public:
    /**
     * Default Constructor.
     */
    LambdaMST3DBy2D() : Super() {}
  };
}// namespace DGtal

#endif // !defined LAMBDAMST3DBy2D_h

#undef LAMBDAMST3DBy2D_RECURSES
#endif // else defined(LAMBDAMST3DBy2D_RECURSES)
