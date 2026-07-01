f=@(x) cos(log(x));
nodes=[0.5 1 1.5 2];
P{1}=@(x) f(nodes(3)).*(x-nodes(4))/(nodes(3)-nodes(4))+f(nodes(4)).*(x-nodes(3))/(nodes(4)-nodes(3));
P{2}=@(x) f(nodes(2)).*(x-nodes(3)).*(x-nodes(4))/((nodes(2)-nodes(3))*(nodes(2)-nodes(4)))+f(nodes(3)).*(x-nodes(2)).*(x-nodes(4))/((nodes(3)-nodes(2))*(nodes(3)-nodes(4)))+f(nodes(4)).*(x-nodes(2)).*(x-nodes(3))/((nodes(4)-nodes(2))*(nodes(4)-nodes(3)));
P{3}=@(x) f(nodes(1)).*(x-nodes(2)).*(x-nodes(3)).*(x-nodes(4))/((nodes(1)-nodes(2)).*(nodes(1)-nodes(3)).*(nodes(1)-nodes(4)))+f(nodes(2)).*(x-nodes(1)).*(x-nodes(3)).*(x-nodes(4))/((nodes(2)-nodes(1))*(nodes(2)-nodes(3))*(nodes(2)-nodes(4)))+f(nodes(3)).*(x-nodes(1)).*(x-nodes(2)).*(x-nodes(4))/((nodes(3)-nodes(1))*(nodes(3)-nodes(2))*(nodes(3)-nodes(4)))+f(nodes(4)).*(x-nodes(1)).*(x-nodes(2)).*(x-nodes(3))/((nodes(4)-nodes(1))*(nodes(4)-nodes(2))*(nodes(4)-nodes(3)));
x=0.4:0.01:2.1;
a=1.75;
figure(1);
plot(x,f(x),'k','LineWidth',2);
hold on;
for i=1:1:3
    plot(x,P{i}(x),'k--','LineWidth',1);
    fprintf("Approximation at the point %.2f: %.8f and Absolute Error: %.8f\n",a,P{i}(a),abs(f(a)-P{i}(a)));
end
plot(nodes,f(nodes),'ro');
hold off;
axis([0.4 2.1 0.5 1.1]);

%%
x=[-0.10;0.00;0.20;0.30];
fx=[17.3000;2.0000;5.1900;1.0000];
coeffs=dividedDiff(x,fx);
points=[0.1;0.4];
for i=1:1:length(points)
    y=dividedDiff_values(points(i),x,coeffs);
    fprintf("Approximated Values at %.2f: %.5f\n",points(i),y);
end

figure(2);
title('Interpolating polynomial');
xlabel('x');
ylabel('Approximations');
plot(x,fx,'ro');
hold on;
data=linspace(-0.2,0.5,100);
values=zeros(1,100);
for i=1:1:100
    values(i)=dividedDiff_values(data(i),x,coeffs);
end
plot(data,values,'k','LineWidth',2);
hold off;

x2=[-0.10;0.00;0.20;0.30;0.05];
fx2=[17.3000;2.0000;5.1900;1.0000;3.1250];
coeffs2=dividedDiff(x2,fx2);
for i=1:1:length(points)
    y2=dividedDiff_values(points(i),x2,coeffs2);
    fprintf("Approximated Values at %.2f: %.5f\n",points(i),y2);
end

figure(3);
title('Interpolating Polynomials');
xlabel('x');
ylabel('Approximated Values');
plot(x2,fx2,'ro');
hold on;
%data=linspace(-0.2,0.5,100);
values2=zeros(1,100);
for i=1:1:100
    values2(i)=dividedDiff_values(data(i),x2,coeffs2);
end
plot(data,values,'k',data,values2,'b','LineWidth',2);
hold off;

function [coeffs]=dividedDiff(x,fx)
table=zeros(length(x));
table(:,1)=fx;
for j=2:1:length(x)
    for i=j:1:length(x)
        table(i,j)=(table(i,j-1)-table(i-1,j-1))/(x(i)-x(i-j+1));
    end
end
coeffs=diag(table);
end

function [y]=dividedDiff_values(p,x,coeffs)
val=coeffs(1);
for k=2:1:length(coeffs)
    pi_k=1;
    for i=1:1:k-1
        pi_k=pi_k*(p-x(i));
    end
    val=val+coeffs(k)*pi_k;
end
y=val;
end

%%
x=[0;0.25;0.5;1;1.25];
a=[0;23.04;47.37;97.45;123.66];
[n,b,c,d]=cubicspline(x,a)
xx=0.75;
i=3;
time=a(i)+b(i)*(xx-x(i))+c(i)*(xx-x(i)).^2+d(i)*(xx-x(i)).^3;
fprintf("Predicted time: %i:%.2f\n",floor(time/60),mod(time,60));

xx=0;
i=1;
speed=3600/(b(i)+2*c(i)*(xx-x(i))+3*d(i)*(xx-x(i)*(xx-x(i))));
fprintf("Predicted speed: %.4f\n",speed);

xx=x(end);
i=n;
speed=3600/(b(i)+2*c(i)*(xx-x(i))+3*d(i)*(xx-x(i))*(xx-x(i)));
fprintf("Predicted speed: %.4f\n",speed);

function [n,b,c,d]=cubicspline(x,a)
n=length(x)-1;
A=zeros(n+1);
h=x(2:end)-x(1:end-1);
A(2:end-1,2:end-1)=diag(2*(h(1:end-1)+h(2:end)))+diag(h(2:end-1),1)+diag(h(2:end-1),-1);
A(1,1)=1;
A(2,1)=h(1);
A(end-1,end)=h(end);
A(end,end)=1;
B=zeros(n+1,1);
B(2:end-1)=3./h(2:end).*(a(3:end)-a(2:end-1))-3./h(1:end-1).*(a(2:end-1)-a(1:end-2));
c=(A\B)';
d=(c(2:end)-c(1:end-1))./(3*h);
b=(a(2:end)-a(1:end-1))./h-h.*(c(2:end)+2*c(1:end-1))/3;
end

