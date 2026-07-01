[e,k]=Error_hilb(10)
plot_error_hilb(3,11);
plot_cond_no_vand(4,17);
cond_no_tridiag(10);
plot_cond_no_tridiag(5,200);

function [e,k]=Error_hilb(n)
A=zeros(n,n);
for i=1:1:n
    for j=1:1:n
        A(i,j)=1/(i+j-1);
    end
end
x_true=ones(n,1);
b=A*x_true;
x_comp=A\b;
x_error=abs(x_true-x_comp);
e=max(x_error);
k=cond(A);
end

function []=plot_error_hilb(n1,n2)
N=n1:n2;
error=zeros(1,n2-n1+1);
condition_number=zeros(1,n2-n1+1);
for i=n1:n2
    [e,k]=Error_hilb(i);
    error(1,i-n1+1)=e;
    condition_number(1,i-n1+1)=k;
end
subplot(1,2,1);
loglog(N,error,'b--.');
xlabel('dimension of matrix');
ylabel('Error');
subplot(1,2,2);
loglog(N,condition_number,'g--*');
xlabel('dimension of matrix');
ylabel('Condition Number K(A)');
title('Vandermonde matrix');
end

function k=cond_no_vand(n)
vand=ones(n+1,n+1);
for i=0:n
    for j=0:n
        vand(i+1,j+1)=(i/n)^j;
    end
end
k=cond(vand);
end

function []=plot_cond_no_vand(n1,n2)
N=n1:n2;
cond=zeros(1,n2-n1+1);
for i=n1:n2
    cond(1,i-n1+1)=cond_no_vand(i);
end
figure;
plot(N,cond,'b--*');
xlabel('Dimension of the matrix');
ylabel('Condition number of the matrix');
title('Vandermonde Matrix');
grid on;
end

function k=cond_no_tridiag(n)
tri=zeros(n,n);
for i=1:1:n
    for j=1:1:n
        if (i==j)
            tri(i,j)=2;
        elseif ((i==(j+1)) || ((i+1)==j))
                tri(i,j)=-1;
        end
    end
end
k=cond(tri);
end

function []=plot_cond_no_tridiag(n1,n2)
N=n1:n2;
cond=zeros(1,n2-n1+1);
for i=n1:n2
    cond(1,i-n1+1)=cond_no_tridiag(i);
end
figure;
plot(N,cond,'b--*');
xlabel('Dimension of the matrix')
ylabel('Condition number of the matrix');
title('Tridiagonal Matrix');
grid on;
end