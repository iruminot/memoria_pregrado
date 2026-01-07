function [ w, x ] = gm_rule_set ( rule, dim_num )

%*****************************************************************************80
%
%% GM_RULE_SET sets a Grundmann-Moeller rule.
%
%  Discussion:
%
%    This is a revised version of the calculation which seeks to compute
%    the value of the weight in a cautious way that avoids intermediate
%    overflow.  Thanks to John Peterson for pointing out the problem on
%    26 June 2008.
%
%    This rule returns weights and abscissas of a Grundmann-Moeller
%    quadrature rule for the DIM_NUM-dimensional unit simplex.
%
%    The dimension POINT_NUM can be determined by calling GM_RULE_SIZE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 June 2008
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Axel Grundmann, Michael Moeller,
%    Invariant Integration Formulas for the N-Simplex
%    by Combinatorial Methods,
%    SIAM Journal on Numerical Analysis,
%    Volume 15, Number 2, April 1978, pages 282-290.
%
%  Parameters:
%
%    Input, integer RULE, the index of the rule.
%    0 <= RULE.
%
%    Input, integer DIM_NUM, the spatial dimension.
%    1 <= DIM_NUM.
%
%    Input, integer POINT_NUM, the number of points in the rule.
%
%    Output, real W(POINT_NUM), the weights.
%
%    Output, real X(DIM_NUM,POINT_NUM), the abscissas.
%
  s = rule;
  d = 2 * s + 1;
  k = 0;
  n = dim_num;
  one_pm = 1;
%point_num = gm_rule_size ( rule, dim_num )
  for i = 0 : s

    weight =  one_pm;

    for j = 1 : max ([ n, d, d + n - i] )

      if ( j <= n )
        weight = weight *  ( j );
      end
      if ( j <= d )
        weight = weight * ( d + n - 2 * i );
      end
      if ( j <= 2 * s )
        weight = weight / 2.0;
      end
      if ( j <= i )
        weight = weight / j;
      end
      if ( j <= d + n - i )
        weight = weight / j;
      end

    end

    one_pm = - one_pm;

    beta_sum = s - i;
    more = 0;
    beta = [];
    h = 0;
    t = 0;

    while ( 1 )

      [ beta, more, h, t ] = comp_next ( beta_sum, dim_num + 1, ...
        beta, more, h, t );

      k = k + 1;

      w(k) = weight;

      x(1:dim_num,k) = ( 2 * beta(2:dim_num+1)' + 1 ) / ( d + n - 2 * i );

      if ( ~more )
        break
      end

    end

  end

  return
end

function point_num = gm_rule_size ( rule, dim_num )

%*****************************************************************************80
%
%% GM_RULE_SIZE determines the size of a Grundmann-Moeller rule.
%
%  Discussion:
%
%    This rule returns the value of POINT_NUM, the number of points associated
%    with a GM rule of given index.
%
%    After calling this rule, the user can use the value of POINT_NUM to
%    allocate space for the weight vector as W(POINT_NUM) and the abscissa
%    vector as X(DIM_NUM,POINT_NUM), and then call GM_RULE_SET.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 July 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Axel Grundmann, Michael Moeller,
%    Invariant Integration Formulas for the N-Simplex
%    by Combinatorial Methods,
%    SIAM Journal on Numerical Analysis,
%    Volume 15, Number 2, April 1978, pages 282-290.
%
%  Parameters:
%
%    Input, integer RULE, the index of the rule.
%    0 <= RULE.
%
%    Input, integer DIM_NUM, the spatial dimension.
%    1 <= DIM_NUM.
%
%    Output, integer POINT_NUM, the number of points in the rule.
%
  arg1 = dim_num + rule + 1;

  point_num = i4_choose ( arg1, rule );

  return
end

function value = i4_choose ( n, k )

%*****************************************************************************80
%
%% I4_CHOOSE computes the binomial coefficient C(N,K).
%
%  Discussion:
%
%    The value is calculated in such a way as to avoid overflow and
%    roundoff.  The calculation is done in integer arithmetic.
%
%    The formula used is:
%
%      C(N,K) = N! / ( K! * (N-K)! )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 June 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    ML Wolfson, HV Wright,
%    Algorithm 160:
%    Combinatorial of M Things Taken N at a Time,
%    Communications of the ACM,
%    Volume 6, Number 4, April 1963, page 161.
%
%  Parameters:
%
%    Input, integer N, K, are the values of N and K.
%
%    Output, integer VALUE, the number of combinations of N
%    things taken K at a time.
%
  mn = min ( k, n - k );

  if ( mn < 0 )

    value = 0;

  elseif ( mn == 0 )

    value = 1;

  else

    mx = max ( k, n - k );
    value = mx + 1;

    for i = 2 : mn
      value = ( value * ( mx + i ) ) / i;
    end

  end

  return
end

function [ a, more, h, t ] = comp_next ( n, k, a, more, h, t )

%*****************************************************************************80
%
%% COMP_NEXT computes the compositions of the integer N into K parts.
%
%  Discussion:
%
%    A composition of the integer N into K parts is an ordered sequence
%    of K nonnegative integers which sum to N.  The compositions (1,2,1)
%    and (1,1,2) are considered to be distinct.
%
%    The routine computes one composition on each call until there are no more.
%    For instance, one composition of 6 into 3 parts is
%    3+2+1, another would be 6+0+0.
%
%    On the first call to this routine, set MORE = FALSE.  The routine
%    will compute the first element in the sequence of compositions, and
%    return it, as well as setting MORE = TRUE.  If more compositions
%    are desired, call again, and again.  Each time, the routine will
%    return with a new composition.
%
%    However, when the LAST composition in the sequence is computed 
%    and returned, the routine will reset MORE to FALSE, signaling that
%    the end of the sequence has been reached.
%
%    This routine originally used a SAVE statement to maintain the
%    variables H and T.  I have decided that it is safer
%    to pass these variables as arguments, even though the user should
%    never alter them.  This allows this routine to safely shuffle
%    between several ongoing calculations.
%
%    There are 28 compositions of 6 into three parts.  This routine will
%    produce those compositions in the following order:
%
%     I         A
%     -     ---------
%     1     6   0   0
%     2     5   1   0
%     3     4   2   0
%     4     3   3   0
%     5     2   4   0
%     6     1   5   0
%     7     0   6   0
%     8     5   0   1
%     9     4   1   1
%    10     3   2   1
%    11     2   3   1
%    12     1   4   1
%    13     0   5   1
%    14     4   0   2
%    15     3   1   2
%    16     2   2   2
%    17     1   3   2
%    18     0   4   2
%    19     3   0   3
%    20     2   1   3
%    21     1   2   3
%    22     0   3   3
%    23     2   0   4
%    24     1   1   4
%    25     0   2   4
%    26     1   0   5
%    27     0   1   5
%    28     0   0   6
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    02 July 2008
%
%  Author:
%
%    Original FORTRAN77 version by Albert Nijenhuis and Herbert Wilf
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Albert Nijenhuis, Herbert Wilf,
%    Combinatorial Algorithms for Computers and Calculators,
%    Second Edition,
%    Academic Press, 1978,
%    ISBN: 0-12-519260-6,
%    LC: QA164.N54.
%
%  Parameters:
%
%    Input, integer N, the integer whose compositions are desired.
%
%    Input, integer K, the number of parts in the composition.
%
%    Input, integer A(K), the previous composition.  On the first call,
%    with MORE = FALSE, set A = [].  Thereafter, A should be the 
%    value of A output from the previous call.
%
%    Input, logical MORE.  The input value of MORE on the first
%    call should be FALSE, which tells the program to initialize.
%    On subsequent calls, MORE should be TRUE, or simply the
%    output value of MORE from the previous call.
%
%    Input, integer H, T, two internal parameters needed for the
%    computation.  The user may need to initialize these before the
%    very first call, but these initial values are not important.
%    The user should not alter these parameters once the computation
%    begins.
%
%    Output, integer A(K), the next composition.
%
%    Output, logical MORE, will be TRUE unless the composition 
%    that is being returned is the final one in the sequence.
%
%    Output, integer H, T, the updated values of the two internal 
%    parameters.
%
  if ( ~more )

    t = n;
    h = 0;
    a(1) = n;
    a(2:k) = 0;

  else
      
    if ( 1 < t )
      h = 0;
    end

    h = h + 1;
    t = a(h);
    a(h) = 0;
    a(1) = t - 1;
    a(h+1) = a(h+1) + 1;

  end

  more = ( a(k) ~= n );

  return
end
