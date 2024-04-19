<?php
/*****************************************************************************/
/* Name : ExponentialRegression.php                                          */
/* Uses : Implementation of several forms for exponential regression.        */
/* Date : 05/06/2013                                                         */
/* Author: Andrew Que <http://www.DrQue.net/>                                */
/* Revisions:                                                                */
/*   1.0 - 05/06/2013 - QUE - Creation.                                      */
/*                                                                           */
/* Implements the following functions:                                       */
/*   ExponentialRegression => a * exp( b * x )                               */
/*   ExponentialRegression2 => a * exp( b * x ) + c                          */
/*  *ExponentialRegression3 => a * exp( b * x + c ) + d                      */
/*   ExponentialDecayRegression => a * ( 1 - exp( b * x ) )                  */
/*   ExponentialDecayRegression2 => a * ( 1 - exp( b * x ) ) + c             */
/*  *ExponentialDecayRegression3 => a * ( 1 - exp( b * x + c ) ) + d         */
/*                                                                           */
/*  * Dysfunctional--Implemented but never converges.                        */
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
 * Various exponential implementation of the Gauss-Newton algorithm.
 *
 * @package GaussNewton
 * @author Andrew Que ({@link http://www.DrQue.net/})
 * @copyright Copyright (c) 2013, Andrew Que
 * @license http://opensource.org/licenses/gpl-license.php GNU Public License
 *
 */

/**
 * f( x ) = a * exp( b * x )
 * Regression for the basic exponential function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentialRegression
 */
class ExponentialRegression extends GaussNewtonRegression
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
    assert( $coefficientIndex < 2 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result = exp( $coefficients[ 1 ] * $x );
        break;
      }

      case 1:
      {
        $result = $coefficients[ 0 ] * $x * exp( $coefficients[ 1 ] * $x );
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
    return 2;
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
    return $coefficients[ 0 ] * exp( $coefficients[ 1 ] * $x );
  }

}

/**
 * f( x ) = a * exp( b * x ) + c
 * Regression for the exponential function with an offset.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentialRegression
 */
class ExponentialRegression2 extends GaussNewtonRegression
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
    assert( $coefficientIndex < 3 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result = exp( $coefficients[ 1 ] * $x );
        break;
      }

      case 1:
      {
        $result = $coefficients[ 0 ] * $x * exp( $coefficients[ 1 ] * $x );
        break;
      }

      case 2:
      {
        $result = 1;
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
    return $coefficients[ 0 ] * exp( $coefficients[ 1 ] * $x ) + $coefficients[ 2 ];
  }

}

/**
 * f( x ) = a * exp( b * x + c ) + d
 * Regression for the exponential function with two offsets.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentialRegression
 */
class ExponentialRegression3 extends GaussNewtonRegression
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
        $result = exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] );

        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] );

        break;
      }

      case 3:
      {
        $result = 1;
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
    return
      $coefficients[ 0 ]
      * exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] )
      + $coefficients[ 3 ];
  }

}

/**
 * f( x ) = a * ( 1 - exp( b * x ) )
 * Regression for the decaying exponential function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentialRegression
 */
class ExponentialDecayRegression extends GaussNewtonRegression
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
    assert( $coefficientIndex < 2 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result = 1 - exp( $coefficients[ 1 ] * $x );
        break;
      }

      case 1:
      {
        $result =
          - $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x );

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
    return 2;
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
    return $coefficients[ 0 ] * ( 1 - exp( $coefficients[ 1 ] * $x ) );
  }

}


/**
 * f( x ) = a * ( 1 - exp( b * x ) ) + c
 * Regression for the decaying exponential function with offset.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentialRegression
 */
class ExponentialDecayRegression2 extends GaussNewtonRegression
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
    assert( $coefficientIndex < 3 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result = 1 - exp( $coefficients[ 1 ] * $x );
        break;
      }

      case 1:
      {
        $result =
          - $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x );

        break;
      }

      case 2:
      {
        $result = 1;

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
    return $coefficients[ 0 ] * ( 1 - exp( $coefficients[ 1 ] * $x ) )
      + $coefficients[ 2 ];
  }

}

/**
 * f( x ) = a * ( 1 - exp( b * x + c ) ) + d
 * Regression for the decaying exponential function with two offsets.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentialRegression
 */
class ExponentialDecayRegression3 extends GaussNewtonRegression
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
        $result = 1 - exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] );
        break;
      }

      case 1:
      {
        $result =
          - $coefficients[ 0 ]
          * $x
          * exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] );

        break;
      }

      case 2:
      {
        $result =
          - $coefficients[ 0 ]
          * exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] );

        break;
      }

      case 3:
      {
        $result = 1;

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
    return
      $coefficients[ 0 ]
      * ( 1 - exp( $coefficients[ 1 ] * $x + $coefficients[ 2 ] ) )
      + $coefficients[ 3 ];
  }

}

?>