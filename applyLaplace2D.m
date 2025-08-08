
function L2 = applyLaplace2D(U, h)
  % 5-točkovni Laplace v 2D (hom. Dirichlet)
  N = size(U,1);
  L2 = zeros(N);
  for i = 1:N
    for j = 1:N
      s = 0;
      if i>1,   s = s + U(i-1,j); end
      if i<N,   s = s + U(i+1,j); end
      if j>1,   s = s + U(i,j-1); end
      if j<N,   s = s + U(i,j+1); end
      L2(i,j) = (s - 4*U(i,j)) / h^2;
    end
  end
end