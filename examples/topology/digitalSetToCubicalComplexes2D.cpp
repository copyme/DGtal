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
 * @file digitalSetToCubicalComplexes2D.cpp
 * @ingroup Examples
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2015/11/03
 *
 * An example file named digitalSetToCubicalComplexes2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include "ConfigExamples.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
// Cellular grid converter
#include "DGtal/topology/DigitalSetToCellularGridConverter.h"
#include "DGtal/topology/helpers/Surfaces.h"
// Shape construction
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"
#include "DGtal/shapes/parametric/Flower2D.h"
// Drawing
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////

int main( int , char** )
{
  trace.beginBlock ( "Example digitalSetToCubicalComplexes2D" );
  trace.beginBlock ( "Generate a 2D shape." );
  typedef Flower2D< Space > MyEuclideanShape;
  MyEuclideanShape shape( RealPoint( 0.0, 0.0 ), 16, 5, 5, M_PI_2/2. );
  
  typedef GaussDigitizer< Space, MyEuclideanShape > MyGaussDigitizer;
  MyGaussDigitizer digShape;
  digShape.attach( shape );
  digShape.init ( shape.getLowerBound(), shape.getUpperBound(), 1.0 );
  Domain domainShape = digShape.getDomain();
  DigitalSet aSet( domainShape );
  Shapes<Domain>::digitalShaper( aSet, digShape );
  
  Board2D board;
  board << SetMode( domainShape.className(), "Paving" ) << domainShape;
  Color dorange ( 255,  136,  0,  220 );
  board << CustomStyle( aSet.className(), new CustomFillColor ( dorange ) );
  board << aSet;
  board.saveEPS ( "digitalSet.eps" );
  board.clear();
  trace.endBlock();
  
  trace.beginBlock ( "Generate a 2D cubical representation." );
  typedef DigitalSetToCellularGridConverter< KSpace, DigitalSet > MyDigitalSetToCellularGridConverter;
  typedef typename MyDigitalSetToCellularGridConverter::Cells Cells;
  typedef typename MyDigitalSetToCellularGridConverter::ConstCellsIterator ConstCellsIterator;
  
  Cells cells0d, cells1d;
  KSpace K;
  K.init (  domainShape.lowerBound(), domainShape.upperBound(), true );
  MyDigitalSetToCellularGridConverter converter ( K );
  converter.init ( aSet.begin(), aSet.end() );
  converter.extractAllCells ( cells0d, 0 );
  converter.extractAllCells ( cells1d, 1 );
  Cell c;
  
  board << SetMode( domainShape.className(), "Paving" ) << domainShape;
  
  board << CustomStyle( Point().className(),
			new CustomColors( Color( 0, 200, 0 ),
					  Color( 100, 255, 100 ) ) );
  for ( typename DigitalSet::ConstIterator it = aSet.begin(); it != aSet.end(); ++it )
    board << *it;
  
  board << CustomStyle( c.className(),
			new CustomColors( Color( 0, 0, 200 ),
					  Color( 100, 100, 255 ) ) );
  
  for ( ConstCellsIterator it = cells0d.begin(); it != cells0d.end(); ++it )
    board << *it;
  board << CustomStyle( c.className(),
			new CustomColors( Color( 200, 0, 0 ),
					  Color( 255, 100, 100 ) ) );
  for ( ConstCellsIterator it = cells1d.begin(); it != cells1d.end(); ++it )
    board << *it;
  
  // Additionally show how to extract boundary.
  SurfelAdjacency<2> sAdj( true );
  vector<vector<SCell> > vectConnectedSCell;
  Surfaces<KSpace>::extractAllConnectedSCell ( vectConnectedSCell,K, sAdj, aSet, false );
  
  for ( unsigned int i = 0; i < vectConnectedSCell.size(); i++ )
    for ( unsigned int j = 0; j < vectConnectedSCell.at(i).size(); j++ )
    {
      board << CustomStyle( vectConnectedSCell.at(i).at(j).className(),
			    new CustomColors( Color( 0, 0, 0 ),
					      Color( 0, 0, 0 ) ) );
      board << vectConnectedSCell.at(i).at(j);
    }
  
  board.saveEPS ( "cubicalComplexes.eps" );
  trace.endBlock();
  trace.endBlock();
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
