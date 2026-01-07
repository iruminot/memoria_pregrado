function n_quad_rule = n_quad_rule(deg,ele_type)
% *********************************************
% ---------------------------------------------------------------
% Description: This function gives the number of points for 
% ************ quadrature rules for exact for polynomials
%             of degree "deg" for triangles (ele_type=2) and edges
%             (ele_type=1).
% ---------------------------------------------------------------
%

% computation of the number of points "n_quad_rule"
% -------------------------------------------------

if (ele_type == 2)  %triangle
    rule = ceil(0.5 * real(deg - 1))	;
    n_quad_rule = gm_rule_size(rule, 2);
elseif (ele_type == 1) % edges
    n_quad_rule  = ceil(0.5 * real(deg + 1))	;	
end


function point_num = gm_rule_size (rule,dim_num)
%**************************************************	
%---------------------------------------------------------------
%	GM_RULE_SIZE determines the size of a Grundmann-Moeller rule.
%
%  	Description:
%    	This rule returns the value of POINT_NUM, the number of points associated
%    	with a GM rule of given index.
%    	After calling this rule, the user can use the value of POINT_NUM to
%    	allocate space for the weight vector as W(POINT_NUM) and the abscissa 
%    	vector as X(DIM_NUM,POINT_NUM), and then call GM_RULE_SET.
%
%  	Parameters:
%    	Input, integer RULE, the index of the rule.
%    	0 <= RULE.
%    	Input, integer DIM_NUM, the spatial dimension.
%    	1 <= DIM_NUM.
%    	Output, integer POINT_NUM, the number of points in the rule.
%---------------------------------------------------------------

arg1 = dim_num + rule + 1;
point_num = binom(arg1, rule);


function value = binom(n, k)
%***************************
%--------------------------------------------------------------
%	Description: The function computes the binomial coefficient 
%			 C(N,K).
%  
%    	The formula used is:
%
%      C(N,K) = N! / (K! * (N-K)!)
%
%    	Parameters:
%    	Input, integer N, K, are the values of N and K.
%       Output, integer binom, the number of combinations of N
%    	things taken K at a time.
%---------------------------------------------------------------

mn = min(k, n - k);
if (mn < 0)
    value = 0;
elseif (mn == 0)
    value = 1;
else
    mx = max(k, n - k);
	value = mx + 1;
	for i = 2 : mn
        value = (value * (mx + i)) / i;
    end
end


