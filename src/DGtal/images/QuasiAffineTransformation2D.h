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
 * @file QuasiAffineTransformation2D.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/09/16
 *
 * Header file for module QuasiAffineTransformation2D.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(QuasiAffineTransformation2D_RECURSES)
#error Recursive header files inclusion detected in QuasiAffineTransformation2D.h
#else // defined(QuasiAffineTransformation2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define QuasiAffineTransformation2D_RECURSES

#if !defined QuasiAffineTransformation2D_h
/** Prevents repeated inclusion of headers. */
#define QuasiAffineTransformation2D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include <functional>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/linalg/SimpleMatrix.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/CImage.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class QuasiAffineTransformation2D
/**
 * Description of class 'QuasiAffineTransformation2D' <p>
 * \brief Aim:
 */

template <typename TSpace>
class QuasiAffineTransformation2D : std::unary_function <typename TSpace::Point, typename TSpace::Point>
{
    ///Checking concepts
    BOOST_CONCEPT_ASSERT(( concepts::CSpace<TSpace> ));
    BOOST_STATIC_ASSERT(( TSpace::dimension == 2 ));
    // ----------------------- Types ------------------------------
public:
    typedef Z2i::Point Point;
    typedef Z2i::RealPoint RealPoint;
    typedef Z2i::Vector Vector;
    typedef SimpleMatrix < int, 2, 2 > Matrix;

    // ----------------------- Interface --------------------------------------
public:
      QuasiAffineTransformation2D ( Matrix matrix, int omega, Z2i::Vector vect )
      {
          m_matrix = matrix;
          m_omega = omega;
          m_vector = vect;
      }
      Point operator() ( const Point & aInput ) const
      {
          RealPoint point ( ( m_matrix * aInput + m_vector ) / m_omega);
          return Point ( std::floor ( point[0] + 0.5f ), std::floor ( point[1] + 0.5f ) );
      }
    // ------------------------- Protected Datas ------------------------------
protected:
    Matrix m_matrix;
    int m_omega;
    Vector m_vector;
};

template <typename TSpace>
class BackwardsQuasiAffineTransformation2D: std::unary_function < typename TSpace::Point, typename TSpace::Point >
{
    ///Checking concepts
    BOOST_CONCEPT_ASSERT(( concepts::CSpace<TSpace> ));
    BOOST_STATIC_ASSERT(( TSpace::dimension == 2 ));
    // ----------------------- Types ------------------------------
public:
    typedef Z2i::Point Point;
    typedef Z2i::RealPoint RealPoint;
    typedef Z2i::Vector Vector;
    typedef std::pair < typename TSpace::Point, typename TSpace::Point > Bounds;
    typedef SimpleMatrix < int, 2, 2 > Matrix;
    // ----------------------- Interface --------------------------------------
public:
    BackwardsQuasiAffineTransformation2D ( const Matrix & matrix, int omega, Z2i::Vector vect )
    {
        int det = matrix.determinant();
        m_matrix(0,0) = matrix(1,1) * omega;
        m_matrix(1,0) = -matrix(1,0) * omega;
        m_matrix(0,1) = -matrix(0,1) * omega;
        m_matrix(1,1) = matrix(0,0) * omega;
        m_vector = ( Matrix() - matrix ) * vect;
        m_omega = std::abs ( det );
    }
    Point operator() ( const Point & aInput ) const
    {
        RealPoint point = ( m_matrix * aInput + m_vector) / m_omega;
        return Point ( std::floor ( point[0] + 0.5f ), std::floor ( point[1] + 0.5f ) );
    }
    // ------------------------- Protected Datas ------------------------------
protected:
    Matrix m_matrix;
    int m_omega;
    Vector m_vector;
};

template <typename TImage>
class QuasiAffineTransformer2D
{
    ///Checking concepts
    BOOST_CONCEPT_ASSERT(( CImage < TImage > ));
    BOOST_STATIC_ASSERT(( TImage::Domain::dimension == 2 ));
    // ----------------------- Types ------------------------------
public:
    typedef Z2i::Point Point;
    typedef Z2i::RealPoint RealPoint;
    typedef Z2i::Vector Vector;
    typedef SimpleMatrix < int, 2, 2 > Matrix;
    // ----------------------- Interface --------------------------------------
public:
    QuasiAffineTransformer2D ( Matrix matrix, int omega, Z2i::Vector vect )
    {
        int det = matrix.determinant();
        m_matrix(0,0) = matrix(1,1) * omega;
        m_matrix(1,0) = -matrix(1,0) * omega;
        m_matrix(0,1) = -matrix(0,1) * omega;
        m_matrix(1,1) = matrix(0,0) * omega;
        m_vector = ( Matrix() - matrix ) * vect;
        m_omega = std::abs ( det );
    }

    // Just calculating starting point and steps
    void BackwardsNN ( const TImage & input, TImage & output )
    {
        Point start = m_matrix * output.domain().lowerBound() + m_vector;
        Vector IncrX = m_matrix * Vector ( 1, - ( output.domain().upperBound()[1] - output.domain().lowerBound()[1] + 1 ) );
        Vector IncrY = m_matrix * Vector ( 0, 1 );

        for ( int i = output.domain().lowerBound()[0]; i <= output.domain().upperBound()[0]; i++ )
        {
            for ( int j = output.domain().lowerBound()[1]; j <= output.domain().upperBound()[1]; j++ )
            {
                Point point ( std::floor ( ( start[0] / m_omega ) + 0.5f ), std::floor ( ( start[1] / m_omega ) + 0.5f ) );
                if ( input.domain().isInside ( point ) )
                    output.setValue ( Point ( i, j ), input ( point ) );
                start += IncrY;
            }
            start += IncrX;
        }
    }

    // ------------------------- Protected Datas ------------------------------
protected:
    Matrix m_matrix;
    int m_omega;
    Vector m_vector;
};

/////////////////////////////////////////////////////////////////////////////
// Template class DomainRigidTransformation2D
/**
     * Description of template functor like class 'DomainRigidTransformation2D' <p>
     * \brief Aim: implements bounds of transformed domain.
     *
     * @tparam TDomain a 2 dimensional domain.
     * @tparam TRigidTransformFunctor a functor which represent two dimensional rigid transformation.
     *
     * @see exampleRigidtransformation2d.cpp
     */
template <typename TDomain, typename TRigidTransformFunctor >
class DomainQuasiAffineTransformation2D :
        std::unary_function < std::pair < typename TDomain::Point, typename TDomain::Point >, TDomain>
{
    ///Checking concepts
    BOOST_STATIC_ASSERT(( TDomain::dimension == 2 ));
    BOOST_CONCEPT_ASSERT(( concepts::CDomain<TDomain> ));

    // ----------------------- Types ------------------------------
public:
    typedef std::pair < typename TDomain::Space::Point, typename TDomain::Space::Point > Bounds;

    // ----------------------- Interface --------------------------------------
public:
    /**
       * Constructor.
       * @param aRigidFunctor  - rigid transformation functor.
       */
    DomainQuasiAffineTransformation2D ( TRigidTransformFunctor & aRigidFunctor ) : transform ( aRigidFunctor ) {}

    /**
       * Operator
       *
       * @return bounds of the transformed domain.
       */
    inline
    Bounds operator()( const TDomain & aInput ) const
    {
        typedef typename TDomain::Point Point;
        Point points[4];
        points[0] = transform ( aInput.lowerBound() );
        points[1] = transform ( aInput.upperBound() );
        points[2] = transform ( Point ( aInput.upperBound()[0], aInput.lowerBound()[1] ) );
        points[3] = transform ( Point ( aInput.lowerBound()[0], aInput.upperBound()[1] ) );

        Point t_min ( INT_MAX, INT_MAX ), t_max ( INT_MIN, INT_MIN );
        for ( unsigned int i = 0; i < 4 ; i++ )
        {
            if ( points[i][0] < t_min[0] )
                t_min[0] = points[i][0];
            if ( points[i][1] < t_min[1] )
                t_min[1] = points[i][1];

            if ( points[i][0] > t_max[0] )
                t_max[0] = points[i][0];
            if ( points[i][1] > t_max[1] )
                t_max[1] = points[i][1];
        }

        Bounds bounds;
        bounds.first = t_min;
        bounds.second = t_max;
        return bounds;
    }

    // ------------------------- Protected Datas ------------------------------
protected:
    TRigidTransformFunctor & transform;
};



} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
//#include "DGtal/images/QuasiAffineTransformation2D.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined QuasiAffineTransformation2D_h

#undef QuasiAffineTransformation2D_RECURSES
#endif // else defined(QuasiAffineTransformation2D_RECURSES)
