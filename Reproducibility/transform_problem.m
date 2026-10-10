function p = transform_problem(original, U, c)
%TRANSFORM_PROBLEM Orthogonal x -> U*x+c for the coupled affine test.
p = original;
p.name = 'coupled_rotated_translated';
p.J = U*original.J*U.';
p.u0 = U*original.u0+c;
p.exact = @(t) bsxfun(@plus,U*original.exact(t),c);
p.forcing = @(t) U*original.forcing(t);
J = p.J; exact = p.exact; F = p.forcing;
p.f = @(t,x) F(t)+J*(x-exact(t));
p.affine_term = @(t) F(t)-J*exact(t);
p.jacobian = @(t,x) J;
end
