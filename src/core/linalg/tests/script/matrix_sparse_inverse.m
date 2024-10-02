function  A_inverse = matrix_sparse_inverse(A, sparsity_pattern)
% Compute sparse inverse matrix of a sparse matrix explicitly.
%
% A               : Matrix to invert
% sparsity_pattern: Sparsity pattern to calculate the sparse inverse of A on
%
% Returns the sparse inverse A^(-1) of the input matrix A.

  n = length(A);
  I = speye(n);
  M = sparse(n,n);

  for j=1:n

    sj=I(:,j);

    col_el=find(sparsity_pattern(j,:));
    A1=A(:,col_el);
    [rows_el,~]=find(A1);
    rows_el=unique(rows_el);
    A1=A1(rows_el,:);
    sj=sj(rows_el,:);

    bj=A1\sj;

    mj = sparse(n,1);
    mj(col_el) = bj;
    A_inverse(:,j) = mj;

  end
end