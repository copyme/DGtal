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
 * @file LambdaMST3D.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/10/06
 *
 * This file is part of the DGtal library.
 */

#if defined(LAMBDAMST3D_RECURSES)
#error Recursive header files inclusion detected in LambdaMST3D.h
#else // defined(LAMBDAMST3D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LAMBDAMST3D_RECURSES

#if !defined LAMBDAMST3D_h
/** Prevents repeated inclusion of headers. */
#define LAMBDAMST3D_h

#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/geometry/curves/CIncrementalSegmentComputer.h"
#include "DGtal/geometry/curves/FunctorsLambdaMST.h"

namespace DGtal {

template < typename TSpace, typename TSegmentation, typename Functor >
class LambdaMST3DEstimator
{
public: 
    ///Checking concepts
    BOOST_CONCEPT_ASSERT(( concepts::CSpace<TSpace> ));
    BOOST_STATIC_ASSERT(( TSpace::dimension == 3 ));
    BOOST_CONCEPT_ASSERT(( CIncrementalSegmentComputer<typename TSegmentation::SegmentComputer> ));
    // ----------------------- Types ------------------------------
public:
    typedef TSegmentation Segmentation;
    typedef typename TSegmentation::SegmentComputer SegmentComputer;
    typedef typename SegmentComputer::ConstIterator ConstIterator;
    typedef typename Functor::Value Value;
    typedef typename TSpace::RealVector RealVector;
    typedef typename TSpace::Point Point;

  // ----------------------- Standard services ------------------------------
public:
  LambdaMST3DEstimator(): myBegin(), myEnd(), dssSegments ( 0 ), myFunctor(), lenFilter ( 0 ), distFilter (-1) {}
    
    /**
     * Initialisation.
     * @param itb, begin iterator
     * @param ite, end iterator
     */
    void init ( const ConstIterator& itb, const ConstIterator& ite )
    {
        myBegin = itb;
        myEnd = ite; 
    }
    
    /**
     * @param segment - DSS segmentation algorithm
     */
    void attach ( const TSegmentation & aSC )
    {
        dssSegments = &aSC;
    }

    void setLengthFilter ( const unsigned int & lenFi )
    {
      lenFilter = lenFi;
    }
    
    //! Works only with  RealVector eval ( const Point & point )
    void setDistanceFilter ( const double & distFi )
    {
      distFilter = distFi;
    }

    /**
     * @param point to calculate A and B for it
     * @return A and B
     */
    RealVector eval ( const Point & point )
    {
      assert ( isValid() );
      typename TSegmentation::SegmentComputerIterator DSS = dssSegments->begin();
      typename TSegmentation::SegmentComputerIterator lastDSS = dssSegments->end();
      Value tangent;
      Value partial;
      Value prev;
      bool found = false;
      bool avg = false;
      
      //!\todo move or remove
      for ( ; DSS != lastDSS; ++DSS )
      {
	if ( !DSS->isInDSS ( point ) && ( EuclideanDistance<Point>( point, *(DSS->begin()) ) <= distFilter || EuclideanDistance <Point>( point, *(DSS->end()-1) ) <= distFilter )  )
	{
	  avg = true;
	  break;
	}
      }
      DSS = dssSegments->begin();
      for ( ; DSS != lastDSS; ++DSS )
      {
	if ( DSS->isInDSS ( point ) || ( EuclideanDistance<Point>( point, *(DSS->begin()) ) <= distFilter || EuclideanDistance <Point>( point, *(DSS->end()-1) ) <= distFilter )  )
	{
	  unsigned int dssLen = std::distance ( DSS.begin(), DSS.end() );
	  if ( dssLen < lenFilter )
	  {
	    std::cerr << "Segment too short: " << dssLen << ", filter: " << lenFilter << std::endl;
	    continue;
	  }
	  
	  prev = partial;
	  found = true;
	  unsigned int pos = std::distance ( DSS.begin(), std::find ( DSS.begin(), DSS.end(), point ) );
	  SegmentComputer comp ( *DSS );
	  partial = myFunctor ( comp, pos, dssLen );
	  if (avg)
	    partial.eccentricity = 1.;
	  if ( prev.vector.cosineSimilarity ( partial.vector ) > M_PI_2 )
	    partial.vector = -partial.vector;
	  partial.vector *= partial.eccentricity;
	  tangent += partial;
	}
      }
      if (!found)
	throw std::runtime_error("Uncovered point!");
      if ( tangent.eccentricity < 1E-20 )
	return tangent.vector;
      else
	return tangent.vector / tangent.eccentricity;
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
      assert ( isValid() );
      std::multimap < Point, Value > outValues;
      typename TSegmentation::SegmentComputerIterator DSS = dssSegments->begin();
      typename TSegmentation::SegmentComputerIterator lastDSS = dssSegments->end();
      
      for(; DSS != lastDSS; ++DSS)
      {
	unsigned int dssLen = std::distance ( DSS.begin(), DSS.end() );
	SegmentComputer comp ( *DSS );
	if ( dssLen > lenFilter )
	  for ( int indexOfPointInDSS = 0; indexOfPointInDSS < dssLen; indexOfPointInDSS++ )
	    outValues.insert ( std::make_pair ( *(DSS.begin() + indexOfPointInDSS), myFunctor ( comp, indexOfPointInDSS, dssLen ) ) );
	else
	  std::cerr << "Segment too short: " << dssLen << ", filter: " << lenFilter << std::endl;
      }
      
      //! \todo simplification -- we need inverse direction when needed, \todo circulator end
      for ( ConstIterator itt = myBegin; itt != myEnd; ++itt )
      {
	typename std::multimap< Point, Value >::const_iterator it  = outValues.lower_bound ( *itt );
	typename std::multimap< Point, Value >::const_iterator it2 = outValues.upper_bound ( *itt );
	static Value prev = it->second;
	Value tangent;
	for (; it != it2; it++ )
	{
	  Value partial = it->second;
	  if (  prev.vector.cosineSimilarity ( partial.vector ) > M_PI_2 )
	  {
	    partial.vector = -partial.vector;
	  }
	  prev = partial;
	  partial.vector *= partial.eccentricity;
	  tangent += partial;
	}
	if ( tangent.eccentricity < 1E-20 )
	  result.push_back ( tangent.vector );
	else
	  result.push_back ( tangent.vector / tangent.eccentricity );
      }
      return;
    }
    
    // ----------------------- Standard services ------------------------------
public:

   /**
   * Checks the validity/consistency of the object.
   * @return 'true' if the object is valid, 'false' otherwise.
   */
    bool isValid() const
    {
        return ( myBegin != myEnd && dssSegments != 0 );
    }
    
    // ------------------------- Internals ------------------------------------
protected:
  
    template <class TPoint>
    double EuclideanDistance ( TPoint const & a, TPoint const & b )
    {
      double dist = 0.;
      for ( unsigned short i = 0; i < TPoint::dimension; i++ )
	dist += ( a[i] - b[i] ) * ( a[i] - b[i] );
      return std::sqrt ( dist );
    }

    // ------------------------- Private Datas --------------------------------
private:
    ConstIterator myBegin;
    ConstIterator myEnd;
    const TSegmentation * dssSegments;
    Functor myFunctor;
    double lenFilter;
    double distFilter;

}; // end of class LambdaTangentFromDSSEstimator 

/**
 * Description of class 'LambdaTangentFromDSS' <p> Aim:
 * \todo -- description
 */
template<typename DSS, typename LambdaFunction>
class TangentFromDSS3DFunctor
{
public:
    // ----------------------- Types ------------------------------
    typedef PointVector<3, double> Vector;
    typedef struct LambdaCharacteristics
    {
        Vector vector;
        double eccentricity = 0.0;
        LambdaCharacteristics & operator += ( const LambdaCharacteristics & ch )
        {
            this->vector += ch.vector;
            this->eccentricity += ch.eccentricity;
            return *this;
        }
    }Value;

    // ----------------------- Interface --------------------------------------
    Value operator() ( const DSS& aDSS, const int & indexOfPointInDSS, const int & dssLen ) const
    {
        Value result;
        typename DSS::Point3d directionZ3;
        Vector direction;
        typename DSS::PointD3d intercept;
        typename DSS::PointD3d thikness;

        aDSS.getParameters ( directionZ3, intercept, thikness );
        direction[0] = directionZ3[0];
        direction[1] = directionZ3[1];
        direction[2] = directionZ3[2];

        result.eccentricity = lambdaFunctor( (double)indexOfPointInDSS / (double)dssLen );

        double norm = direction.norm();
        if ( norm > 1E-20 )
            direction /= norm;
        result.vector = direction;
        return result;
    }
private:
    // ------------------------- Private Datas --------------------------------
    LambdaFunction lambdaFunctor;
};

//-------------------------------------------------------------------------------------------

// Template class LambdaMST3D
/**
* \brief Aim: Simplify creation of Lambda MST tangent estimator.
*
*/
template < typename DSSSegmentationComputer, typename lambdaFunction = functors::Lambda64Function>
class LambdaMST3D:
        public LambdaMST3DEstimator<Z3i::Space, DSSSegmentationComputer,
        TangentFromDSS3DFunctor< typename DSSSegmentationComputer::SegmentComputer, lambdaFunction> >
{
    typedef 
    LambdaMST3DEstimator<Z3i::Space, DSSSegmentationComputer,
    TangentFromDSS3DFunctor< typename DSSSegmentationComputer::SegmentComputer, lambdaFunction> > Super;
    
public: 
    /**
   * Default Constructor.
   */
    LambdaMST3D() : Super() {}
};
}// namespace DGtal

#endif // !defined LAMBDAMST3D_h

#undef LAMBDAMST3D_RECURSES
#endif // else defined(LAMBDAMST3D_RECURSES)
