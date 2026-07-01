
y=Newton(1,1e-8,50);
function [y]=f(x)
y=x*x-2*exp(-x)*x+exp(-2*x);
end
function [y]=df(x)
y=2*x-2*exp(-x)+2*x*exp(-x)-2*exp(-2*x);
end

function [y]= Newton(x0,eps,maxItr)
x_prev=x0;
err=Inf;
currItr=0;
while(err>eps && currItr<=maxItr)
    x_new=x_prev-f(x_prev)/df(x_prev);
    err=abs((x_new-x_prev)/x_new);
    currItr=currItr+1;
    x_prev=x_new;
end
fprintf("Got the answer after %d iterations and the answer is %.8f\n",currItr,x_new);
y=x_new;
end
