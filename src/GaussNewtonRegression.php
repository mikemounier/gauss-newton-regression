<?php
/*****************************************************************************/
/* Name : GaussNewtonRegression.php                                          */
/* Uses : Class for preforming Gauss-Newton regression.                      */
/* Date : 05/05/2013                                                         */
/* Author: Andrew Que <http://www.DrQue.net/>                                */
/* Revisions:                                                                */
/*   1.0 - 05/05/2013 - QUE - Creation.                                      */
/*                                                                           */
/* Credits:                                                                  */
/*   Prof. C. Balaji of the Indian Institute of Technology Madras.           */
/*   http://mech.iitm.ac.in/Faculty/CB/home.php                              */
/*                                                                           */
/* References:                                                               */
/*   (1)...http://en.wikipedia.org/wiki/Matrix_multiplication                */
/*   (2)...http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm       */
/*   (3)...http://en.wikipedia.org/wiki/Gaussian_elimination                 */
/*   (4)...http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant      */
/*   (5)...http://en.wikipedia.org/wiki/Coefficient_of_determination         */
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

/**
 * Gauss-Newton non-linear regression.
 *
 * Used for calculating least-square regression coefficients.  Useful for
 * linear and non-linear curve fitting.
 *
 * @package GaussNewton
 * @author Andrew Que ({@link http://www.DrQue.net/})
 * @copyright Copyright (c) 2013, Andrew Que
 * @license http://opensource.org/licenses/gpl-license.php GNU Public License
 *
 */

/**
 * Gauss-Newton non-linear regression.
 *   This class is useful for computing the coefficients for some non-linear
 * equation using the Gauss-Newton algorithm.
 *
 * Algorithm (see 2):
 *     b_n+1 = b_n + (J^T J)^-1 J r( b_n )
 *   Where b is the current guess for the coefficients, n is the current guess
 * J is the Jacobian, and r is the residual (error).
 *
 * This is an abstract class and can not be used by itself.  It requires the
 * the function and function partial differentials to be define.
 *
 * @package GaussNewtonRegression
 * @abstract
 */
abstract class GaussNewtonRegression
{
  // Implementation notes:
  //   All matrices in this class are assumed to be a two-dimensional array
  // indexed by row and then column.  Matricies are only used internally, and
  // all public input/output is done using one-dimensional arrays.
  //   Matrix inversion is not implemented because the multiplying the a matrix
  // by the an other matrix's inverse is the same as solving a system of
  // equations.  Solving this system directly saves the multiplication step.
  //

  /**
   * Compute partial-differential of x.
   *
   * Compute the partial-differential with respect to a given coefficient for
   * a given value of x.
   *     p( x ) = d f( x ) / d c_n
   * Example:
   * <code>
   *   If f( x ) = c_0 e^(c_1 x ), then:
   *     d f( x ) / d c_0 = e^( c_1 x )
   *     d f( x ) / d c_1 = c_0 x e^( c_1 x )
   *   Each function must be defined.
   * </code>
   * @param float $x Value of x to supply to partial-differential function.
   * @param integer $coefficientIndex Which coefficient to be used in the
   *      partial-differential function.
   * @param array $coefficients Values of the coefficients to be used.
   * @return float p( x ) for the supplied input.
   */
  abstract protected function partialDifferential
  (
    $x,
    $coefficientIndex,
    $coefficients
  );

  /**
   * Number of coefficients.
   *
   * @return integer The number of coefficients for the equation being used.
   * @abstract
   */
  abstract protected function getNumberOfCoefficients();

  /**
   * Evaluate function at x.
   *
   * Return f( x ) for a given set of coefficients.
   *
   * @param float $x Real value given to equation.
   * @param array $coefficients Coefficients used in calculation.
   * @return float Value result of equation.
   * @abstract
   */
  abstract public function getFunction( $x, $coefficients );

  /**
   * Return the transpose of a matrix.
   *
   * @param array $matrix Any matrix.
   * @return array A matrix representing transpose.
   */
  private function transpose( $matrix )
  {
    $result = array();
    for ( $row = 0; $row < count( $matrix[ 0 ] ); ++$row )
    {
      $result[ $row ] = array();
      for ( $column = 0; $column < count( $matrix ); ++$column )
        $result[ $row ][ $column ] = $matrix[ $column ][ $row ];
    }

    return $result;
  }

  /**
   * Multiply two matrices (see 1).
   *   [ A ] = [ B ][ C ]
   *
   * @param array $matrix - Matrix [ B ].
   * @param array $multiplicand - Matrix [ C ].
   * @return array  Matrix [ C ].
   */
  private function multiply( $matrix, $multiplicand )
  {
    $result = array();

    for ( $row = 0; $row < count( $matrix ); ++$row )
    {
      $result[ $row ] = array();
      for ( $column = 0; $column < count( $multiplicand[ 0 ] ); ++$column )
      {
        $result[ $row ][ $column ] = 0;

        for ( $index = 0; $index < count( $matrix[ $row ] ); ++$index )
          $result[ $row ][ $column ] +=
            $matrix[ $row ][ $index ] * $multiplicand[ $index ][ $column ];
      }
    }

    return $result;
  }

  /**
   * Get a matrix of how much error is between the calculated value, and the
   * true value.
   *
   * @param float $x Array of x-coordinates.
   * @param float $y Array of known y-coordinates that correspond to $x.
   * @param array $coefficients Current coefficients to check for error.
   * @return array Matrix of error at each point.
   */
  private function getErrorMatrix( $x, $y, $coefficients )
  {
    assert( count( $x ) == count( $y ) );

    $result = array();
    for ( $row = 0; $row < count( $x ); ++$row )
      $result[ $row ] =
        array( $y[ $row ] - $this->getFunction( $x[ $row ], $coefficients ) );

    return $result;
  }

  /**
   * Solve a system of equations in matrix form.
   *
   * Done by concatenating answer matrix to the left, doing Gaussian
   * elimination (see 3), and returning the last column.
   *    [ M ][ C ] = [ A ] ==> [ C ] = [ M ]^1 [ A ]
   *
   * @param array $matrix - Square matrix [ M ].
   * @param array $answers - Single column matrix [ A ].
   * @output array Single column matrix [ C ].
   */
  private function solve( $matrix, $answers )
  {
    $degree = count( $matrix );
    $order = array();

    //-----------------------------------
    // Add the answers to the matrix.
    //-----------------------------------
    $isDone = array();
    for ( $row = 0; $row < $degree; ++$row )
    {
      $matrix[ $row ][ $degree ] = $answers[ $row ][ 0 ];
      $isDone[ $row ] = false;
    }

    //-----------------------------------
    // We now row-reduce the matrix to obtain the coefficients for the
    // polynomial.
    //-----------------------------------

    // This loop will result in an upper-triangle matrix with the
    // diagonals all 1--the first part of row-reduction--using 2
    // elementary row operations: multiplying a row by a scalar, and
    // subtracting a row by a multiple of an other row.
    // NOTE: This loop can be done out-of-order.  That is, the first
    // row may not begin with the first term.  Order is tracked in the
    // "Order" array.
    $order = array();
    for ( $column = 0; $column < $degree; ++$column )
    {
      // Find a row to work with.
      // A row that has a term in this column, and hasn't yet been
      // reduced.
      $activeRow = 0;
      while ( ( ( 0 == $matrix[ $activeRow ][ $column ] )
             || ( $isDone[ $activeRow ] ) )
           && ( $activeRow < $degree ) )
      {
        ++$activeRow;
      }

      // Do we have a term in this row?
      if ( $activeRow < $degree )
      {
        // Remeber the order.
        $order[ $column ] = $activeRow;

        // Normilize row--results in the first term being 1.
        $firstTerm = $matrix[ $activeRow ][ $column ];
        for ( $subColumn = $column; $subColumn <= $degree; ++$subColumn )
          $matrix[ $activeRow ][ $subColumn ] /= $firstTerm;

        // This row is finished.
        $isDone[ $activeRow ] = true;

        // Subtract the active row from all rows that are not finished.
        for ( $row = 0; $row < $degree; ++$row )
          if ( ( ! $isDone[ $row ] )
            && ( 0 != $matrix[ $row ][ $column ] ) )
          {
             // Get first term in row.
             $firstTerm = $matrix[ $row ][ $column ];
             for ( $subColumn = $column; $subColumn <= $degree; ++$subColumn )
               $matrix[ $row ][ $subColumn ] -=
                 $firstTerm * $matrix[ $activeRow ][ $subColumn ];
          }
      }
    }

    // Reset done.
    for ( $row = 0; $row < $degree; ++$row )
      $isDone[ $row ] = false;

    // Storage for the resulting coefficients.
    $coefficients = array();

    // Back-substitution.
    // This will solve the matrix completely, resulting in the identity
    // matrix in the x-locations, and the coefficients in the last column.
    //   | 1  0  0 ... 0  c0 |
    //   | 0  1  0 ... 0  c1 |
    //   | .  .  .     .   . |
    //   | .  .  .     .   . |
    //   | 0  0  0 ... 1  cn |
    for ( $column = ($degree - 1); $column >= 0; --$column )
    {
      // The active row is based on order.
      $activeRow = $order[ $column ];

      // The active row is now finished.
      $isDone[ $activeRow ] = true;

      // For all rows not finished...
      for ( $row = 0; $row < $degree; ++$row )
        if ( ! $isDone[ $row ] )
        {
          $firstTerm = $matrix[ $row ][ $column ];

          // Back substation.
          for ( $subColumn = $column; $subColumn <= $degree; ++$subColumn )
            $matrix[ $row ][ $subColumn ] -=
              $firstTerm * $matrix[ $activeRow ][ $subColumn ];
        }

      // Save this coefficients for the return.
      $coefficients[ $column ] = array( $matrix[ $activeRow ][ $degree ] );
    }

    return $coefficients;
  }

  /**
   * Refine coefficients one around.
   *
   * Run an iteration of the Gauss-Newton algorithm.  Run this as many times
   * as needed to refine coefficients to their best values.
   *
   * @param array $x Array of x data points.
   * @param array $y Array of y data points that correspond to x.
   * @param array $coefficients Initial or previous guess at coefficients.
   * @return array An array that is the new guess at coefficients.
   */
  public function refineCoefficients( $x, $y, $coefficients )
  {
    // Number of coefficients.
    $numberOfCoefficients = $this->getNumberOfCoefficients();

    // Compute Jacobian matrix (see 4).
    // The Jacobian matrix consists of the partial-differentials.
    //  [ d f( x_0 ) / d c_0   d f( x_0 ) / d c_1  ...  d f( x_0 ) / d c_m ]
    //  [ d f( x_1 ) / d c_0   d f( x_1 ) / d c_1  ...  d f( x_1 ) / d c_m ]
    //  [       .                  .              .           .            ]
    //  [       .                  .                .         .            ]
    //  [       .                  .                  .       .            ]
    //  [ d f( x_n ) / d c_0   d f( x_n ) / d c_1  ...  d f( x_n ) / d c_m ]
    //
    //   Where f() is the non-linear function, d f() / d c is the partial-
    // differential with respect to c.  n is the number of rows in x, and m
    // is the number of coefficients in the equation.
    $jacobianMatrix = array();
    for ( $row = 0; $row < count( $x ); ++$row )
    {
      $jacobianMatrix[ $row ] = array();
      for ( $column = 0; $column < $numberOfCoefficients; ++$column )
        $jacobianMatrix[ $row ][ $column ] =
          $this->partialDifferential( $x[ $row ], $column, $coefficients );
    }

    // Compute the transpose of Jacobian matrix.
    $jacobianTranspose = $this->transpose( $jacobianMatrix );

    // Compute J^T * J.
    // Where J is the Jacobian matrix, and J^T is the transpose.
    $jacobianIntermediate =
      $this->multiply( $jacobianTranspose, $jacobianMatrix );

    // Get the residual error R.  This matrix defines how "good" the fit is.
    $errorMatrix = $this->getErrorMatrix( $x, $y, $coefficients );

    // Compute J^T * R.
    $errorIntermediate = $this->multiply( $jacobianTranspose, $errorMatrix );

    // Compute the amount of change needed to each coefficient.
    $coefficientsDelta =
      $this->solve( $jacobianIntermediate, $errorIntermediate );

    // Adjust coefficient by change needed.
    for ( $index = 0; $index < $this->getNumberOfCoefficients(); ++$index )
      $coefficients[ $index ] += $coefficientsDelta[ $index ][ 0 ];

    // Return new guess at coefficients.
    return $coefficients;
  }

  /**
   * Compute R-Squared.
   *
   * Compute Coefficient of determination (R squared) for a set of
   * coefficients.  This is useful for determining how good a fit the
   * coefficients are.  (See 5)
   *
   * @param array $x Array of x data points.
   * @param array $y Array of y data points that correspond to x.
   * @param array $coefficients Current coefficients to compare against.
   * @return float Real value between 0 and 1.
   */
  public function getR_Squared( $x, $y, $coefficients )
  {
    $errorMatrix = $this->getErrorMatrix( $x, $y, $coefficients );

    $average = array_sum( $y ) / count( $y );

    $totalSumOfSquares = 0;
    $residualSumOfSquares = 0;
    foreach ( $errorMatrix as $index => $error )
    {
      $totalSumOfSquares += pow( $y[ $index ] - $average, 2 );
      $residualSumOfSquares += $error[ 0 ] * $error[ 0 ];
    }

    return 1 - ( $residualSumOfSquares / $totalSumOfSquares );
  }

} // GaussNewtonRegression

?>