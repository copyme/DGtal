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
#include <vector>
#include <map>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/geometry/curves/CIncrementalSegmentComputer.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/curves/LambdaMST2D.h"

namespace DGtal {

  template < typename TSegmentationXY, typename TSegmentationXZ, typename TSegmentationYZ, typename Iterator, typename Functor >
class LambdaMST3DBy2DEstimator
{
public:
    ///Checking concepts
    BOOST_CONCEPT_ASSERT(( CIncrementalSegmentComputer<typename TSegmentationXY::SegmentComputer> ));
    BOOST_CONCEPT_ASSERT(( CIncrementalSegmentComputer<typename TSegmentationXZ::SegmentComputer> ));
    BOOST_CONCEPT_ASSERT(( CIncrementalSegmentComputer<typename TSegmentationYZ::SegmentComputer> ));
    // ----------------------- Types ------------------------------
public:
    typedef TSegmentationXY SegmentationXY;
    typedef TSegmentationXZ SegmentationXZ;
    typedef TSegmentationYZ SegmentationYZ;
    typedef typename TSegmentationXY::SegmentComputer SegmentComputerXY;
    typedef typename TSegmentationXZ::SegmentComputer SegmentComputerXZ;
    typedef typename TSegmentationYZ::SegmentComputer SegmentComputerYZ;
    typedef PointVector<3, double> RealVector3D;
    typedef PointVector<3, int> Point3D;
    typedef PointVector<2, int> Point2D;
    typedef typename Functor::MAIN_AXIS MAIN_AXIS;
    typedef functors::Projector < SpaceND < 2, int > > Projector2d;

  // ----------------------- Standard services ------------------------------
public:
  LambdaMST3DBy2DEstimator(): myBegin(), myEnd(), dssSegmentsXY ( 0 ), dssSegmentsXZ ( 0 ), dssSegmentsYZ ( 0 ) 
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
    void init ( const Iterator& itb, const Iterator& ite )
    {
        myBegin = itb;
        myEnd = ite; 
    }
    
    /**
     * @param segment - DSS segmentation algorithm
     */
    void attach ( const TSegmentationXY & aSCXY, const TSegmentationXZ & aSCXZ, const TSegmentationYZ & aSCYZ )
    {
        dssSegmentsXY = &aSCXY;
        dssSegmentsXZ = &aSCXZ;
        dssSegmentsYZ = &aSCYZ;
	
	lmst64XY.attach ( *dssSegmentsXY );
	lmst64XZ.attach ( *dssSegmentsXZ );
	lmst64YZ.attach ( *dssSegmentsYZ );
    }

    /**
     * @param point to calculate A and B for it
     * @return A and B
     */
    RealVector3D eval ( const Point3D & point )
    {
      assert ( isValid() );
      MAIN_AXIS axis = detectMainAxis ( point );
      RealVector3D tangent;
      typename LambdaMST2D<SegmentationXY>::RealVector v0;
      typename LambdaMST2D<SegmentationXY>::RealVector v1;
      lmst64XY.eval ( myProjXY ( point ) );
      if ( axis == MAIN_AXIS::XY )
      {
	v0 = lmst64XY.eval ( myProjXY ( point ) );
	v1 = lmst64XZ.eval ( myProjXZ ( point ) );
      }
      else if ( axis == MAIN_AXIS::XZ )
      {
	v0 = lmst64XY.eval ( myProjXY ( point ) );
	v1 = lmst64YZ.eval ( myProjYZ ( point ) );
      }
      else
      {
	v0 = lmst64XZ.eval ( myProjXZ ( point ) );
	v1 = lmst64YZ.eval ( myProjYZ ( point ) );
      }
      return myFunctor ( axis, v0, v1 );
    }

    /**
     * @param result output iterator on the estimated quantity
     *
     * @return the estimated quantity
     * from itb till ite (excluded)
     */
    template <typename ResultType>
    void eval ( std::vector<ResultType> & result )
    {
      assert ( myBegin != myEnd && isValid() );
      for( Iterator it = myBegin; it != myEnd; ++it )
      {
	RealVector3D tan = eval ( *it );
	if ( result.size() > 0 && tan.cosineSimilarity ( result.back() ) > M_PI_2 )
	  tan -= tan;
	result.push_baCk ( tan );
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
        return ( dssSegmentsXY && dssSegmentsXZ && dssSegmentsYZ );
    }
    
    // ------------------------- Internals ------------------------------------
protected:
  inline
  MAIN_AXIS detectMainAxis ( const Point3D & p )
  {
    unsigned int rankXY = 0, rankXZ = 0, rankYZ = 0;
    for ( typename SegmentationXY::SegmentComputerIterator it = dssSegmentsXY->begin(); it != dssSegmentsXY->end(); ++it )
      if ( it->isInDSS ( myProjXY ( p ) ) )
	rankXY += std::distance ( it->begin(), it->end() );
    for ( typename SegmentationXZ::SegmentComputerIterator it = dssSegmentsXZ->begin(); it != dssSegmentsXZ->end(); ++it )
      if ( it->isInDSS ( myProjXZ ( p ) ) )
	rankXZ += std::distance ( it->begin(), it->end() );
    for ( typename SegmentationYZ::SegmentComputerIterator it = dssSegmentsYZ->begin(); it != dssSegmentsYZ->end(); ++it )
      if ( it->isInDSS ( myProjYZ ( p ) ) )
	rankYZ += std::distance ( it->begin(), it->end() );
      
      if ( rankXY > rankYZ && rankXZ > rankYZ )
	return MAIN_AXIS::XY;
      else if ( rankXY > rankXZ && rankYZ > rankXZ )
	return MAIN_AXIS::XZ;
      else
	return MAIN_AXIS::YZ;
  }

    // ------------------------- Private Datas --------------------------------
private:
    Iterator myBegin;
    Iterator myEnd;
    const TSegmentationXY * dssSegmentsXY;
    const TSegmentationXZ * dssSegmentsXZ;
    const TSegmentationYZ * dssSegmentsYZ;
    LambdaMST2D < SegmentationXY > lmst64XY;
    LambdaMST2D < SegmentationXZ > lmst64XZ;
    LambdaMST2D < SegmentationYZ > lmst64YZ;
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
     enum MAIN_AXIS {XY = 0, XZ = 1, YZ = 2};
     
     // ----------------------- Interface --------------------------------------
     Vector3D operator() ( MAIN_AXIS mainAxis, const Vector2D & v0, const Vector2D & v1 ) const
     {
       Vector3D tangent;
       if ( mainAxis == XY )
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
       else if ( mainAxis == XZ )
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
typename DSSSegmentationComputerXY,
typename DSSSegmentationComputerXZ, 
typename DSSSegmentationComputerYZ, 
typename Iterator >
class LambdaMST3DBy2D:
public LambdaMST3DBy2DEstimator < DSSSegmentationComputerXY, DSSSegmentationComputerXZ, DSSSegmentationComputerYZ, Iterator,
  TangentFromDSS3DBy2DFunctor >
  {
    typedef LambdaMST3DBy2DEstimator < DSSSegmentationComputerXY, DSSSegmentationComputerXZ, DSSSegmentationComputerYZ, Iterator,
    TangentFromDSS3DBy2DFunctor > Super;
    
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
