<?php
/*****************************************************************************/
/* Name : ExponentiationRegression.php                                       */
/* Uses : Implementation of several power regressions functions.             */
/* Date : 05/06/2013                                                         */
/* Author: Andrew Que <http://www.DrQue.net/>                                */
/* Revisions:                                                                */
/*   0.1 - 05/06/2013 - QUE - Creation.                                      */
/*   1.0 - 05/23/2013 - QUE - Several additional methods added.              */
/*                                                                           */
/* Implements the following functions:                                       */
/*   f( x ) = a * x^b                                                        */
/*   f( x ) = a * x^b + c                                                    */
/*   f( x ) = a * n^( b * x )                                                */
/*   f( x ) = a * n^( b * x ) + c                                            */
/*   f( x ) = a * b^( x )                                                    */
/*  *f( x ) = a * b^( c * x )                                                */
/*  *f( x ) = a * b^( c * x + d )                                            */
/*  *f( x ) = a * b^( c * x + d ) + g                                        */
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
 * Various implementations of exponentiation regression.
 *
 * @package GaussNewton
 * @author Andrew Que ({@link http://www.DrQue.net/})
 * @copyright Copyright (c) 2013, Andrew Que
 * @license http://opensource.org/licenses/gpl-license.php GNU Public License
 *
 */

/**
 * f( x ) = a * x^b
 * Regression for the basic exponentiation function.
 *
 * Note that this function can also be used for inverse functions in the
 * form f( x ) = a / x^b be having a negative value for b.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class ExponentiationRegression extends GaussNewtonRegression
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
        $result =
          pow
          (
            $x,
            $coefficients[ 1 ]
          );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * pow
            (
              $x,
              $coefficients[ 1 ]
            )
          * log( $x );

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
    return $coefficients[ 0 ] * pow( $x, $coefficients[ 1 ] );
  }

}

/**
 * f( x ) = a * x^b + c
 * Regression for the basic exponentiation function.
 *
 * Note that this function can also be used for inverse functions in the
 * form f( x ) = a / x^b + c be having a negative value for b.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class ExponentiationRegression2 extends GaussNewtonRegression
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
        $result = pow( $x, $coefficients[ 1 ] );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * pow( $x, $coefficients[ 1 ] )
          * log( $x );

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
    return
      $coefficients[ 0 ]
      * pow( $x, $coefficients[ 1 ] )
      + $coefficients[ 2 ];
  }

}

/**
 * f( x ) = a * n^( b * x )
 * Regression for the basic exponentiation function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class FixedPowerRegression extends GaussNewtonRegression
{
  /**
   * The base power to use.
   */
  private $n;

  /**
   * Constructor.
   *
   * @param float $n Base power.
   */
  public function __construct( $n )
  {
    $this->n = $n;
  }

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
        $result = pow( $this->n, $coefficients[ 1 ] * $x );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * $this->n
          * pow( $this->n, $coefficients[ 1 ] * $x )
          * log( $this->n );

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
    return $coefficients[ 0 ] * pow( $this->n, $coefficients[ 1 ] * $x  );
  }
}

/**
 * f( x ) = a * n^( b * x ) + c
 * Regression for the basic exponentiation function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class FixedPowerRegression2 extends GaussNewtonRegression
{
  /**
   * The base power to use.
   */
  private $n;

  /**
   * Constructor.
   *
   * @param float $n Base power.
   */
  public function __construct( $n )
  {
    $this->n = $n;
  }

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
        $result = pow( $this->n, $coefficients[ 1 ] * $x );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * $x
          * pow( $this->n, $coefficients[ 1 ] * $x )
          * log( $this->n );

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
    return
      $coefficients[ 0 ]
      * pow( $this->n, $coefficients[ 1 ] * $x )
      + $coefficients[ 2 ] ;
  }

}

/**
 * f( x ) = a * b^( x )
 * Regression for the basic exponentiation function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class PowerRegression extends GaussNewtonRegression
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
        $result = pow( $coefficients[ 1 ], $x );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * pow( $coefficients[ 1 ], $x - 1 )
          * $x;

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
    return
      $coefficients[ 0 ] * pow( $coefficients[ 1 ], $x );
  }

}

/**
 * f( x ) = a * b^( c * x )
 * Regression for the basic exponentiation function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class PowerRegression2 extends GaussNewtonRegression
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
        $result =
          pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x
          );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * pow
            (
              $coefficients[ 1 ],
              $coefficients[ 2 ] * $x - 1
            )
          * $coefficients[ 2 ] * $x;

        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * pow
            (
              $coefficients[ 1 ],
              $coefficients[ 2 ] * $x
            )
          * $x
          * log( $coefficients[ 1 ] );

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
    return
      $coefficients[ 0 ] * pow( $coefficients[ 1 ], $coefficients[ 2 ] * $x );
  }

}

/**
 * f( x ) = a * b^( c * x + d )
 * Regression for the basic exponentiation function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class PowerRegression3 extends GaussNewtonRegression
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
          pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x + $coefficients[ 3 ]
          );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * pow
            (
              $coefficients[ 1 ],
              $coefficients[ 2 ] * $x + $coefficients[ 3 ] - 1
            )
          * ( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );

        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * pow
            (
              $coefficients[ 1 ],
              $coefficients[ 2 ] * $x + $coefficients[ 3 ]
            )
          * $x
          * log( $coefficients[ 1 ] );

        break;
      }

      case 3:
      {
        $result =
          $coefficients[ 0 ]
          * pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x + $coefficients[ 3 ]
          )
          * log( $coefficients[ 1 ] );

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
      * pow( $coefficients[ 1 ], $coefficients[ 2 ] * $x + $coefficients[ 3 ] );
  }

}

/**
 * f( x ) = a * b^( c * x + d ) + g
 * Regression for the basic exponentiation function.
 *
 * @package GaussNewtonRegression
 * @subpackage ExponentiationRegression
 */
class PowerRegression4 extends GaussNewtonRegression
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
    assert( $coefficientIndex < 5 );

    $result = 0;
    switch ( $coefficientIndex )
    {
      case 0:
      {
        $result =
          pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x + $coefficients[ 3 ]
          );
        break;
      }

      case 1:
      {
        $result =
          $coefficients[ 0 ]
          * pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x + $coefficients[ 3 ] - 1
          )
          * ( $coefficients[ 2 ] * $x + $coefficients[ 3 ] );

        break;
      }

      case 2:
      {
        $result =
          $coefficients[ 0 ]
          * pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x + $coefficients[ 3 ]
          )
          * $x
          * log( $coefficients[ 1 ] );

        break;
      }

      case 3:
      {
        $result =
          $coefficients[ 0 ]
          * pow
          (
            $coefficients[ 1 ],
            $coefficients[ 2 ] * $x + $coefficients[ 3 ]
          )
          * log( $coefficients[ 1 ] );

        break;
      }

      case 4:
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
    return 5;
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
      * pow( $coefficients[ 1 ], $coefficients[ 2 ] * $x + $coefficients[ 3 ] )
      + $coefficients[ 4 ];
  }

}

?>