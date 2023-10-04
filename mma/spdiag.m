%SPDIAGS Sparse matrix formed from diagonals.
%   SPDIAGS, which generalizes the function "diag", deals with three
%   matrices, in various combinations, as both input and output.
%
%   [B,d] = SPDIAGS(A) extracts all nonzero diagonals from the m-by-n
%   matrix A.  B is a min(m,n)-by-p matrix whose columns are the p
%   nonzero diagonals of A.  d is a vector of length p whose integer
%   components specify the diagonals in A.
%
%   B = SPDIAGS(A,d) extracts the diagonals specified by d.
%   A = SPDIAGS(B,d,A) replaces the diagonals of A specified by d with
%       the columns of B.  The output is sparse.
%   A = SPDIAGS(B,d,m,n) creates an m-by-n sparse matrix from the
%       columns of B and places them along the diagonals specified by d.
%
%   Roughly, A, B and d are related by
%       for k = 1:p
%           B(:,k) = diag(A,d(k))
%       end
%
%   Example: These commands generate a sparse tridiagonal representation
%   of the classic second difference operator on n points.
%       e = ones(n,1);
%       A = spdiags([e -2*e e], -1:1, n, n)
%
%   Some elements of B, corresponding to positions "outside" of A, are
%   not actually used.  They are not referenced when B is an input and
%   are set to zero when B is an output.  If a column of B is longer than
%   the diagonal it's representing, elements of super-diagonals of A
%   correspond to the lower part of the column of B, while elements of
%   sub-diagonals of A correspond to the upper part of the column of B.
%
%   Example: This uses the top of the first column of B for the second
%   sub-diagonal and the bottom of the third column of B for the first
%   super-diagonal.
%       B = repmat((1:n)',1,3);
%       S = spdiags(B,[-2 0 1],n,n);
%
%   See also DIAG, SPEYE.

%   Rob Schreiber
%   Copyright 1984-2019 The MathWorks, Inc.

function [res1,res2] = spdiag(arg1,arg2,arg3,arg4)

if nargin <= 2
    % Extract diagonals
    A = arg1;
    if nargin == 1
        % Find all nonzero diagonals
        [i,j] = find(A);
        % Compute d = unique(d) without extra function call
        d = sort(j-i);
        d = d(diff([-inf; d(:)])~=0);
        d = d(:);
    else
        % Diagonals are specified
        d = arg2(:);
    end
    [m,n] = size(A);
    p = length(d);
    B = zeros(min(m,n),p,class(A));
    for k = 1:p
        if m >= n
            i = max(1,1+d(k)):min(n,m+d(k));
        else
            i = max(1,1-d(k)):min(m,n-d(k));
        end
        B(i,k) = diagk(A,d(k));
    end
    res1 = B;
    res2 = d;
end

if nargin >= 3
    B = arg1;
    d = arg2(:);
    p = length(d);
    if nargin == 3 % Replace specified diagonals
        A = arg3;
    else           % Create new matrix with specified diagonals
        A = sparse(arg3, arg4);
    end
    [m,n] = size(A);
    
    % Check size of matrix B (should be min(m,n)-by-p)
    % For backwards compatibility, only error if the code would
    % previously have errored out in the indexing expression.
    maxIndexRows = max(max(1,1-d), min(m,n-d)) + (m>=n)*d;
    maxIndexRows(max(1,1-d) > min(m,n-d)) = 0;
    if any(maxIndexRows > size(B, 1)) || p > size(B, 2)
        if nargin == 3
            error(message('MATLAB:spdiags:InvalidSizeBThreeInput'));
        else
            error(message('MATLAB:spdiags:InvalidSizeBFourInput'));
        end
    end
    
    % Compute indices and values of sparse matrix with given diagonals
    if nnz(A) == 0 && size(B, 1) == min(m, n) && size(B, 2) == p && p > 1
        if m < n
            % Compute transpose of A, then transpose before returning.
            d = -d;
        end
        
        % Sort d in descending order and reorder B accordingly:
        if issorted(d, 'ascend')
            d = flip(d);
            B = flip(B, 2);
        elseif ~issorted(d, 'descend')
            [d, ind] = sort(d, 'descend');
            B = B(:, ind);
        end
        
        if m >= n
            res1 = makeSparsePresorted(B, d, m, n);
        else
            % Compute transpose of A, then transpose before returning.
            res1 = makeSparsePresorted(B, d, n, m);
            res1 = res1.';
        end
    else
        % Insert diagonals into existing matrix A. This algorithm
        % has fewer requirements on B than the one above.
        res1 = makeSparseGeneral(B, d, A);
    end
    
    if islogical(A) || islogical(B)
        res1 = (res1~=0);
    end
end
end


function A = makeSparsePresorted(B, d, m, n)
% Compute inputs to the sparse function in such a way
% that they will all be sorted, for performance.
% The input d must be in descending order.

% Compute arrays of lower and upper bounds for each column
% The nonzeros of column A(:, ii) are contained in B(lb(ii):ub(ii), ii).
% If A(:, ii) is all zero, then either lb(ii) or ub(ii) is NaN.
% Compact definition:
%   lb(ii) = find(ii > d(:), 1, 'first');
%   ub(ii) = find(ii <= d(:) + m, 1, 'last');
lb = nan(1, n);
if ~isempty(d)
    ii = n;
    for lowerBound=1:length(d)
        while ii >= 1 && ii > d(lowerBound)
            lb(ii) = lowerBound;
            ii = ii-1;
        end
    end
end

ub = nan(1, n);
if ~isempty(d)
    ii = 1;
    for upperBound=length(d):-1:1
        while ii <= n && ii <= d(upperBound) + m
            ub(ii) = upperBound;
            ii = ii+1;
        end
    end
end

% Number of nonzeros in each column
blocklen = ub-lb+1;

% If either ub(ii) or lb(ii) is NaN, then A(:, ii) is all zero
% and should be skipped.
blocklen(isnan(blocklen)) = 0;

% Column index
J = repelem(1:n, blocklen);
J = J(:);

% Row index
I = zeros(length(J), 1);
offset = 0;
for ii=1:n
    for jj=1:blocklen(ii)
        I(offset+jj) = ii - d(lb(ii) + jj - 1);
    end
    offset = offset + blocklen(ii);
end

% Values vector
V = zeros(length(J), 1);
offset = 0;
for ii=1:n
    for jj=1:blocklen(ii)
        V(offset+jj) = B(ii, lb(ii) + jj - 1);
    end
    offset = offset + blocklen(ii);
end

A = sparse(I, J, V, m, n);
end


function A = makeSparseGeneral(B, d, A)
% Construct sparse matrix by inserting
% diagonals into existing matrix A.

[m, n] = size(A);
p = length(d);

% Precompute number of nonzeros (parts of each diagonal that overlap with
% the matrix) and allocate inputs I, J and V for sparse
nz = sum(max(0, min(m, n-d) - max(1, 1-d) + 1));
I = zeros(nz, 1);
J = zeros(nz, 1);
V = zeros(nz, 1);

% Fill in the diagonals
offset = 1;
for k = 1:p
    % Append new d(k)-th diagonal to compact form
    for i=max(1, 1-d(k)):min(m, n-d(k))
        I(offset) = i;
        J(offset) = i + d(k);
        if m >= n
            V(offset) = B(i + d(k), k);
        else
            V(offset) = B(i, k);
        end
        offset = offset + 1;
    end
end

if nnz(A) > 0
    % Process A in compact form
    [Iold,Jold,Vold] = find(A);
    
    % Delete current d(k)-th diagonal, k=1,...,p
    i = any((Jold(:) - Iold(:)) == d', 2);
    Iold(i) = [];
    Jold(i) = [];
    Vold(i) = [];
    
    % Combine new diagonals and non-diagonal entries of original matrix
    I = [I(:); Iold(:)];
    J = [J(:); Jold(:)];
    V = [V(:); Vold(:)];
end

A = sparse(I, J, V, m, n);
end


function D = diagk(X,k)
% DIAGK  K-th matrix diagonal.
% DIAGK(X,k) is the k-th diagonal of X, even if X is a vector.
if ~isvector(X)
    D = diag(X,k);
    D = D(:);  %Ensure column vector is returned for empty X.
else
    if ~isempty(X) && 0 <= k && 1+k <= size(X,2)
        D = X(1+k);
    elseif ~isempty(X) && k < 0 && 1-k <= size(X,1)
        D = X(1-k);
    else
        D = zeros(0,1,'like',X);
    end
end
end

