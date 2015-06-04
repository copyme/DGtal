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
 * @file LambdaMST2D.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/10/03
 *
 * This file is part of the DGtal library.
 */

#if defined(LAMBDAMST2D_RECURSES)
#error Recursive header files inclusion detected in LambdaMST2D.h
#else // defined(LAMBDAMST2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LAMBDAMST2D_RECURSES

#if !defined LAMBDAMST2D_h
/** Prevents repeated inclusion of headers. */
#define LAMBDAMST2D_h

#include <algorithm>
#include <cmath>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/geometry/curves/CIncrementalSegmentComputer.h"
#include "DGtal/geometry/curves/FunctorsLambdaMST.h"

/**
 * Aim: Implement Lambda MST tangent estimator.
 */

namespace DGtal {
template < typename TSpace, typename TSegmentation, typename Functor >
class LambdaMST2DEstimator
{
    ///Checking concepts
    BOOST_CONCEPT_ASSERT(( concepts::CSpace<TSpace> ));
    BOOST_STATIC_ASSERT(( TSpace::dimension == 2 ));
    BOOST_CONCEPT_ASSERT(( CIncrementalSegmentComputer<typename TSegmentation::SegmentComputer> ));
    // ----------------------- Types ------------------------------
public: 
    typedef TSegmentation Segmentation;
    typedef typename TSegmentation::SegmentComputer SegmentComputer;
    typedef typename SegmentComputer::ConstIterator ConstIterator;
    typedef typename Functor::Value Value;
    typedef typename TSpace::RealVector RealVector;
    typedef typename TSpace::Point Point;
    
    // ----------------------- Interface --------------------------------------
public:
    LambdaMST2DEstimator() : dssSegments ( 0 ) { }
    
    /**
     * Initialisation.
     * @param itb, begin iterator
     * @param ite, end iterator
     */
    void init ( const ConstIterator & itb, const ConstIterator & ite )
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
    
    /**
     * @param point to calculate A and B for it
     * @return A and B
     */
    RealVector eval ( const Point & point )
    {
      assert ( dssSegments != 0 );
      typename TSegmentation::SegmentComputerIterator DSS = dssSegments->begin();
      typename TSegmentation::SegmentComputerIterator lastDSS = dssSegments->end();
      Value tangent;
      bool found = false;
      for ( ; DSS != lastDSS; ++DSS )
      {
	if ( DSS->isInDSS ( point ) )
	{
	  found = true;
	  unsigned int pos = std::distance ( DSS.begin(), std::find ( DSS.begin(), DSS.end(), point ) );
	  unsigned int dssLen = std::distance ( DSS.begin(), DSS.end() );
	  SegmentComputer comp ( *DSS );
	  tangent += myFunctor ( comp, pos, dssLen );
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
     * More efficient way to compute tangent for all points in a curve.
     * @param Containter to store results.
     */
    void eval ( std::vector <RealVector> & result )
    {
      assert ( isValid() );
      std::multimap < Point, Value > outValues;
      typename TSegmentation::SegmentComputerIterator DSS = dssSegments->begin();
      typename TSegmentation::SegmentComputerIterator lastDSS = dssSegments->end();
      for ( ; DSS != lastDSS; ++DSS )
      {
	unsigned int dssLen = std::distance ( DSS.begin(), DSS.end() );
	SegmentComputer comp ( *DSS );
	for ( int indexOfPointInDSS = 0; indexOfPointInDSS < dssLen; indexOfPointInDSS++ )
	  outValues.insert ( std::make_pair ( *(DSS.begin() + indexOfPointInDSS), myFunctor ( comp, indexOfPointInDSS, dssLen ) ) );
      }
      for ( ConstIterator itt = myBegin; itt != myEnd; ++itt )
      {
	typename std::multimap< Point, Value >::const_iterator it  = outValues.lower_bound ( *itt );
	typename std::multimap< Point, Value >::const_iterator it2 = outValues.upper_bound ( *itt );
	Value tangent;
	for (; it != it2; ++it )
	  tangent += it->second;
	if ( tangent.eccentricity < 1E-20 )
	  result.push_back ( tangent.vector );
	else
	  result.push_back ( tangent.vector / tangent.eccentricity );
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
        return ( myBegin != myEnd && dssSegments != 0 );
    }

    // ------------------------- Private Datas --------------------------------
private:

    ConstIterator myBegin;
    ConstIterator myEnd;
    Functor myFunctor;
    const TSegmentation * dssSegments;

}; // end of class LambdaTangentFromDSSEstimator 

/**
 * Description of class 'LambdaTangentFromDSS' <p> Aim:
 * \todo -- description
 */
template<typename DSS, typename LambdaFunction>
class TangentFromDSS2DFunctor
{
    // ----------------------- Types ------------------------------
    typedef PointVector<2, double> Vector;
public:
    typedef struct LambdaCharacteristics
    {
        Vector vector;
        double eccentricity = 0.0;
        LambdaCharacteristics & operator+= ( const LambdaCharacteristics & ch )
        {
            this->vector += ch.vector;
            this->eccentricity += ch.eccentricity;
            return *this;
        }
    }Value;

    // ----------------------- Interface --------------------------------------
public:
    Value operator() ( const DSS& aDSS, const int & indexOfPointInDSS, const int & dssLen ) const
    {
        Value result;
        double norm = std::sqrt ( aDSS.a() * aDSS.a() + aDSS.b() * aDSS.b() );
        result.eccentricity = lambdaFunctor( (double)indexOfPointInDSS / (double)dssLen );
        if ( norm > 1E-20 )
        {
            result.vector[0] = result.eccentricity * aDSS.a () / norm;
            result.vector[1] = result.eccentricity * aDSS.b () / norm;
        }
        else
        {
            result.vector[0] = result.eccentricity * aDSS.a ();
            result.vector[1] = result.eccentricity * aDSS.b ();
        }
        return result;
    }
private:
    // ------------------------- Private Datas --------------------------------
    LambdaFunction lambdaFunctor;
};

//-------------------------------------------------------------------------------------------

// Template class LambdaMST2D
/**
* \brief Aim: Simplify creation of Lambda MST tangent estimator.
*
*/
template < typename DSSSegmentationComputer, typename LambdaFunction = functors::Lambda64Function >
class LambdaMST2D: 
        public LambdaMST2DEstimator < Z2i::Space, DSSSegmentationComputer,
        TangentFromDSS2DFunctor < typename DSSSegmentationComputer::SegmentComputer, LambdaFunction > >
{
    typedef LambdaMST2DEstimator < Z2i::Space, DSSSegmentationComputer,
    TangentFromDSS2DFunctor < typename DSSSegmentationComputer::SegmentComputer, LambdaFunction> > Super;
    
public: 
    /**
    * Default Constructor.
    */
    LambdaMST2D (): Super() {}
};
}// namespace DGtal

#endif // !defined LAMBDAMST2D_h

#undef LAMBDAMST2D_RECURSES
#endif // else defined(LAMBDAMST2D_RECURSES)
