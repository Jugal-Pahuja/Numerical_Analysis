y1=roots(6,5,-4)
y2=roots(6e154,5e154,-4e154)
y3=roots(0,1,1)
y4=roots(1,-1e5,1)
y5=roots(1,-4,3.999999)
y6=roots(1e-155,1e155,1e155)
y7=roots(1,1,1)

function y=roots(a,b,c)
s=max([abs(a),abs(b),abs(c)]);
if s==0
    error("All the values are zero; equation is indeterminate");
end
a=a/s;
b=b/s;
c=c/s;
if a==0
    root=-c/b;
    y=[0,root,root];
    return;
end
det=(b^2)-4*a*c;
if det<0
    % Case 1: Complex roots
    x1=(-b)/(2*a);
    x2=sqrt(-det)/(2*a);
    y=[1,x1,x2];
else
    % Case 2:Real roots
    if b<0
        x1=(-b+sqrt(det))/(2*a);
        x2=(2*c)/(-b+sqrt(det));
    else
        x1=(-b-sqrt(det))/(2*a);
        x2=(2*c)/(-b-sqrt(det));
    end
    y=[1,x1,x2];
end
end

function y=fin_diff_tan(x,h)
y=(tan(x+h)-tan(x))/h;
end

function e=error_fd_tan(x,h)
% absolute error
e=abs((sec(x)^2)-fin_diff_tan(x,h));
end

function y=cent_diff_tan(x,h)
y=(tan(x+h)-tan(x-h))/(2*h);
end

function e=error_cd_tan(x,h)
e=abs((sec(x)^2)-cent_diff_tan(x,h));
end

function []=plot_error(x,K)
k=0:K;
plot(10.^(-k),error_fd_tan(x,10.^(-k)));
end

