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

/**
 * @file testLambdaMST2D.cpp
 * @ingroup Tests
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/10/03
 *
 * Functions for testing class LambdaMST2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <vector>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/LambdaMST2D.h"
#include "DGtal/curves/ParametricCurveDigitizer3D.h"
#include "DGtal/curves/parametric/EllipticHelix.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#ifdef __GNUC__
   #ifndef NDEBUG
      #include <fenv.h>
   #endif
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class LambdaMST2D.
///////////////////////////////////////////////////////////////////////////////

class testLambdaMST3D
{
    typedef EllipticHelix < Z3i::Space > MyHelix;
    typedef DGtal::Z2i::Point Point;
    typedef std::vector < Point > Range;
    typedef Range::const_iterator ConstIterator;
    typedef ParametricCurveDigitizer3D < MyHelix >  Digitizer;
    typedef ParametricCurveDigitizer3D < MyHelix >::DigitalCurve MyDigitalCurve;
    typedef std::vector<Z2i::Point> MyProjection;
    typedef ArithmeticalDSSComputer < ConstIterator, int, 8 > SegmentComputer;
    typedef SaturatedSegmentation<SegmentComputer> Segmentation;
    typedef functors::Projector< SpaceND<2,int> > Projector2d;
private:
    Digitizer digitizer;
    MyHelix helix;
    MyDigitalCurve digitalCurve;
    MyProjection oXY;
public:
    testLambdaMST3D () : helix ( 10, 20, 0 ) {
        digitizer.attach ( helix );
        digitizer.init ( digitalCurve );
	double maxt = M_PI * 2.;
        digitizer.digitize ( 0, maxt, 0.001 );

        std::vector<DGtal::Dimension> v1;
        v1.push_back(0);
        v1.push_back(1);
        Projector2d myProjXY;
        myProjXY.init ( v1.begin(), v1.end() );

        for ( typename MyDigitalCurve::const_iterator it = digitalCurve.begin(); it != digitalCurve.end(); ++it )
        {
            if ( oXY.size() == 0 )
                oXY.push_back ( myProjXY ( *it ) );
            else if ( oXY.back() != myProjXY ( *it ) )
                oXY.push_back ( myProjXY ( *it ) );
        }
    }
    bool lambda64ByPoint ()
    {
        Segmentation segmenter ( oXY.cbegin(), oXY.cend(), SegmentComputer() );
        LambdaMST2D < Segmentation > lmst64;
        lmst64.attach ( segmenter );
        lmst64.init ( oXY.begin(), oXY.end() );
        for ( unsigned int i = 0; i < oXY.size(); i++ )
        {
            lmst64.eval ( oXY[i] );
        }
        return true;
    }
    bool lambda64()
    {
        Segmentation segmenter ( oXY.cbegin(), oXY.cend(), SegmentComputer() );
        LambdaMST2D < Segmentation > lmst64;
        lmst64.attach ( segmenter );
        lmst64.init ( oXY.begin(), oXY.end() );
        std::vector < LambdaMST2D < Segmentation >::RealVector > tangent;
        lmst64.eval ( tangent );
        return true;
    }
    bool lambda64Both()
    {
        Segmentation segmenter ( oXY.cbegin(), oXY.cend(), SegmentComputer() );
        LambdaMST2D < Segmentation > lmst64;
        lmst64.attach ( segmenter );
        lmst64.init ( oXY.begin(), oXY.end() );
        std::vector < LambdaMST2D < Segmentation >::RealVector > tangent;
        lmst64.eval ( tangent );
        for ( unsigned int i = 0; i < oXY.size(); i++ )
        {
            if ( lmst64.eval ( oXY[i] ) != tangent[i] )
                return false;
        }
        return true;
    }
};


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int , char**  )
{
#ifdef __GNUC__
   #ifndef NDEBUG
    fetestexcept ( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW );
   #endif
#endif
    bool res = true;
    testLambdaMST3D testLMST;
    trace.beginBlock ( "Testing LambdaMST2D" );
        trace.beginBlock ( "Testing point only calculation" );
          res &= testLMST.lambda64ByPoint();
        trace.endBlock();
        trace.beginBlock ( "Testing calculation for whole curve" );
           res &= testLMST.lambda64();
        trace.endBlock();
        trace.beginBlock ( "Testing values obtined from the both methods." );
           res &= testLMST.lambda64Both();
        trace.endBlock();
    trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
    trace.endBlock();
    return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
