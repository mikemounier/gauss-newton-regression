<?php
/*****************************************************************************/
/* Name : SineRegression.php                                                 */
/* Uses : Implementation of several forms sine regression.                   */
/* Date : 05/14/2013                                                         */
/* Author: Andrew Que <http://www.DrQue.net/>                                */
/* Revisions:                                                                */
/*   1.0 - 05/14/2013 - QUE - Creation.                                      */
/*                                                                           */
/* Implements the following functions:                                       */
/*   f( x ) = a * exp( b * x ) * sin( c * x + d )                            */
/*   f( x ) = a * exp( b * x ) * ( cos( c * x ) + sin( c * x ) )             */
/*   f( x ) = a * exp( b * x ) * ( cos( c * x + d ) + sin( c * x + d ) )     */
/*                                                                           */
/* License:                                                                  */
/*   This program is free software: you can redistribute it and/or modify    */
/*   it under the terms of the GNU General Public License as published by    */
/*   the Free Software Foundation, either version 3 of the License, or       */
/*   (at your option) any later version.                                     */
/*                                                                           */
/*   This program is distributed in the hope that it will be useful,         */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*   GNU General Public License for more details.                            */
/*                                                                           */
/*   You should have received a copy of the GNU General Public License       */
/*   along with this program.  If not, see <http://www.gnu.org/licenses/>.   */
/*                                                                           */
/*                     (C) Copyright 2013 by Andrew Que                      */
/*                           http://www.DrQue.net/                           */
/*****************************************************************************/

require_once( "GaussNewtonRegression.php" );

/**
 * Various exponentiation (power-of) implementation of the Gauss-Newton
 * algorithm.
 *
 * @package GaussNewton
 * @author Andrew Que ({@link http://www.DrQue.net/})
 * @copyright Copyright (c) 2013, Andrew Que
 * @license http://opensource.org/licenses/gpl-license.php GNU Public License
 *
 */

/**
 * f( x ) = a * exp( b * x ) * sin( c * x + d )
 * Regression for the dampened sine function.
 *
 * @package GaussNewtonRegression
 * @subpackage SineRegression
 */
class DampedSineRegression extends GaussNewtonRegression
{
  /**
   * Partial differential.
   *
   * @param float $x Value of x to supply to partial-differential function.
   * @param integer $coefficientIndex Which coefficient to be used in the
   *      partial-differential function.
   * @param array $coefficients Values of the coefficients to be used.
   * @return float p( x ) for the supplied input.
   */
  protected function partialDifferential( $x, $coefficientIndex, $coefficients )
  {
    assert( $coefficientIndex < 4 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result =
          exp( $coefficients[ 1 ] * $x )
          * sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x )
          * sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );
        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x )
          * cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );
        break;
      }

      case 3:
      {
        $result =
          $coefficients[ 0 ]
          * exp( $coefficients[ 1 ] * $x )
          * cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );
        break;
      }

      default:
      {
        assert( false );
        break;
      }
    }

    return $result;
  }


  /**
   * Number of coefficients.
   *
   * @return integer The number of coefficients for the equation being used.
   */
  protected function getNumberOfCoefficients()
  {
    return 4;
  }

  /**
   * Evaluate function at x.
   *
   * Return f( x ) for a given set of coefficients.
   *
   * @param float $x Real value given to equation.
   * @param array $coefficients Coefficients used in calculation.
   * @return float Value result of equation.
   */
  public function getFunction( $x, $coefficients )
  {
    return $coefficients[ 0 ]
            * exp( $coefficients[ 1 ] * $x )
            * sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );
  }

}

/**
 * f( x ) = a * exp( b * x ) * ( cos( c * x ) + sin( c * x ) )
 * Regression for the full dampened sine function.
 *
 * @package GaussNewtonRegression
 * @subpackage SineRegression
 */
class DampedSineRegression2 extends GaussNewtonRegression
{
  /**
   * Partial differential.
   *
   * @param float $x Value of x to supply to partial-differential function.
   * @param integer $coefficientIndex Which coefficient to be used in the
   *      partial-differential function.
   * @param array $coefficients Values of the coefficients to be used.
   * @return float p( x ) for the supplied input.
   */
  protected function partialDifferential( $x, $coefficientIndex, $coefficients )
  {
    assert( $coefficientIndex < 4 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result =
          exp( $coefficients[ 1 ] * $x )
          * ( cos( $coefficients[ 2 ] * $x )
            + sin( $coefficients[ 2 ] * $x ) );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x )
          * ( cos( $coefficients[ 2 ] * $x )
            + sin( $coefficients[ 2 ] * $x ) );
        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * exp( $coefficients[ 1 ] * $x )
          * ( $x * cos( $coefficients[ 2 ] * $x )
            - $x * sin( $coefficients[ 2 ] * $x ) );
        break;
      }

      default:
      {
        assert( false );
        break;
      }
    }

    return $result;
  }

  /**
   * Number of coefficients.
   *
   * @return integer The number of coefficients for the equation being used.
   */
  protected function getNumberOfCoefficients()
  {
    return 3;
  }

  /**
   * Evaluate function at x.
   *
   * Return f( x ) for a given set of coefficients.
   *
   * @param float $x Real value given to equation.
   * @param array $coefficients Coefficients used in calculation.
   * @return float Value result of equation.
   */
  public function getFunction( $x, $coefficients )
  {
    return $coefficients[ 0 ]
            * exp( $coefficients[ 1 ] * $x )
            * ( cos( $coefficients[ 2 ] * $x )
              + sin( $coefficients[ 2 ] * $x ) );
  }

}

/**
 * f( x ) = a * exp( b * x ) * ( cos( c * x + d ) + sin( c * x + d ) )
 * Regression for the full dampened sine function with phase.
 *
 * @package GaussNewtonRegression
 * @subpackage SineRegression
 */
class DampedSineRegression3 extends GaussNewtonRegression
{
  /**
   * Partial differential.
   *
   * @param float $x Value of x to supply to partial-differential function.
   * @param integer $coefficientIndex Which coefficient to be used in the
   *      partial-differential function.
   * @param array $coefficients Values of the coefficients to be used.
   * @return float p( x ) for the supplied input.
   */
  protected function partialDifferential( $x, $coefficientIndex, $coefficients )
  {
    assert( $coefficientIndex < 4 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result =
          exp( $coefficients[ 1 ] * $x )
          * ( cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] )
            + sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] ) );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x )
          * ( cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] )
            + sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] ) );
        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x )
          * ( cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] )
            - sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] ) );
        break;
      }

      case 3:
      {
        $result =
          $coefficients[ 0 ]
          * exp( $coefficients[ 1 ] * $x )
          * ( cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] )
            - sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] ) );
        break;
      }

      default:
      {
        assert( false );
        break;
      }
    }

    return $result;
  }

  /**
   * Number of coefficients.
   *
   * @return integer The number of coefficients for the equation being used.
   */
  protected function getNumberOfCoefficients()
  {
    return 4;
  }

  /**
   * Evaluate function at x.
   *
   * Return f( x ) for a given set of coefficients.
   *
   * @param float $x Real value given to equation.
   * @param array $coefficients Coefficients used in calculation.
   * @return float Value result of equation.
   */
  public function getFunction( $x, $coefficients )
  {
    return $coefficients[ 0 ]
            * exp( $coefficients[ 1 ] * $x )
            * ( cos( $coefficients[ 2 ] * $x + $coefficients[ 3 ] )
              + sin( $coefficients[ 2 ] * $x + $coefficients[ 3 ] ) );
  }

}
?>