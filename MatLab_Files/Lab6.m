Bisection_Method(0,1,1e-8);
Bisection_Method(0,1,1e-12);

x=[0]*20;
for i=-100:1:100
    x(i+101)=i;
end
y1=x;
y2=2*sin(x);
plot(x,y1,x,y2);

function y=func(x)
y=x-(2^(-x));
end

function []=Bisection_Method(a0,b0,eps)
max_itr=30;
err=(b0-a0)/2;
itr=1;
a=a0;
b=b0;
while(err>=eps && itr<=max_itr)
    x=(a+b)/2;
    if (func(x)==0)
        fprintf("Found the root %d",x);
        break;
    end
    if (sign(func(a))*sign(func(x))<0)
        b=x;
    end
    if (sign(func(x))*sign(func(b))<0)
        a=x;
    end
    err=(b-a)/2;
    itr=itr+1;
end
fprintf("Tolerance: %d, Approxiamation: %.12f, Iterations: %d\n",eps,x,itr);
end


function []=Bisection_Method_part2(a0,b0,eps)
max_itr=30;
err=(b0-a0)/2;
itr=1;
a=a0;
b=b0;
while(err>=eps && itr<=max_itr)
    x=(a+b)/2;
    if (func(x)==0)
        fprintf("Found the root %d",x);
        break;
    end
    if (sign(func(a))*sign(func(x))<0)
        b=x;
    end
    if (sign(func(x))*sign(func(b))<0)
        a=x;
    end
    err=(b-a)/2;
    itr=itr+1;
end
fprintf("Tolerance: %d, Approxiamation: %.12f, Iterations: %d\n",eps,x,itr);
end


