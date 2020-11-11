syms x y z

eqn = cos(y)^2+sin(y)^2-z^2+z*z-5*y+6*x;

tic
simplify(eqn,'Steps',3);
toc