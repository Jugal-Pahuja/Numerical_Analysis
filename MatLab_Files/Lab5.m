b=[6;25;-11;15];
A=[10,-1,2,0; -1,11,-1,3; 2,-1,10,-1; 0,3,-1,8];
x0=[0;0;0;0];
[y,itr]=jacobi(x0,A,b,1e-4,50);
fprintf("The soltuion of Jacobi method is [%.4f, %.4f, %.4f, %.4f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),itr-1);
[y,itr]=gaussSeidel(x0,A,b,1e-4,50);
fprintf("The soltuion of Gauss Seidel method is [%.4f, %.4f, %.4f, %.4f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),itr-1);
[J,GS]=matrix(A);
fprintf("Spectral radius of Jacobi matrix is: %.6f\n",J);
fprintf("Spectral radius of Gauss-Seidel matrix is: %.6f\n",GS);

b=[0;5;0;6;2;-6];
x0=[0;0;0;0;0;0];
eps=1e-8;
max_itr=50;
A=[4,-1,0,0,0,0; -1,4,-1,0,0,0; 0,-1,4,0,0,0; 0,0,0,4,-1,0; 0,0,0,-1,4,-1; 0,0,0,0,-1,4];
[y,itr]=jacobi(x0,A,b,1e-8,50);
fprintf("The soltuion of Jacobi method is [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),y(5),y(6),itr-1);
[y,itr]=gaussSeidel(x0,A,b,1e-8,50);
fprintf("The soltuion of Gauss Seidel method is [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),y(5),y(6),itr-1);

[y,itr]=sor(1.01,A,b,eps,x0,max_itr);
fprintf("The soltuion of SOR method is [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),y(5),y(6),itr-1);
[y,itr]=sor(1.05,A,b,eps,x0,max_itr);
fprintf("The soltuion of SOR method is [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),y(5),y(6),itr-1);
[y,itr]=sor(1.10,A,b,eps,x0,max_itr);
fprintf("The soltuion of SOR method is [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),y(5),y(6),itr-1);
[y,itr]=sor(1.15,A,b,eps,x0,max_itr);
fprintf("The soltuion of SOR method is [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f] and the number of iterations being %d\n",y(1),y(2),y(3),y(4),y(5),y(6),itr-1);

[J,GS]=matrix(A);
w_opt=2/(1+sqrt(1-J^2));
fprintf("The optimal value of omega is %.4f\n",w_opt);

function [y,itr]=jacobi(x0,A,b,eps,max_itr)
err=inf;
itr=1;
x=x0;
D=diag(diag(A));
P=A-D;
while itr<=max_itr && err>=eps
    x=D\(b-P*x);
    err=norm(b-A*x,2)/norm(b,2);
    itr=itr+1;
end
y=x;
end

function [y,itr]=gaussSeidel(x0,A,b,eps,max_itr)
x=x0;
itr=1;
err=inf;
N=tril(A);
U=A-N;
while itr<=max_itr && err>=eps
    x=N\(b-U*x);
    err=norm(b-A*x,2)/norm(b,2);
    itr=itr+1;
end
y=x;
end

function [J,GS]=matrix(A)
D=diag(diag(A));
L=tril(A,-1);
U=triu(A,1);
M_J=-D\(L+U);
J=max(abs(eig(M_J)));
M_GS=-(L+D)\U;
GS=max(abs(eig(M_GS)));
end

function [y,itr]=sor(w,A,b,eps,x0,max_itr)
itr=1;
err=inf;
x=x0;
D=diag(diag(A));
L=tril(A,-1);
U=triu(A,1);
while itr<=max_itr && err>=eps
    x=(D+w*L)\((D*(1-w)-w*U)*x+w*b);
    err=norm(b-A*x,2)/norm(b,2);
    itr=itr+1;
end
y=x;
end
