function [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcnApax(obj, t, x, T, controller, Grad1, Grad2)
%EventFcnApax dy=0

fun = x(7);

if ~Grad1 
    grad_x   = 0;
    grad_xi  = 0;
    grad_eps = 0;
    grad_t   = 0;
    grad_T   = 0;
else
    grad_x  = [zeros(1,6) 1 0 0 0];
    grad_xi = 0*zeros(1,controller.nXi);

    grad_eps = 0;

    grad_t   = 0;
    grad_T   = 0;
end
end

