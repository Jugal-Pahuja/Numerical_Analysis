A1=[2 0 1 -1;6 3 2 -1;4 3 -2 3;-2 -6 2 -14];
b1=[6;15;3;12];
%x1=Gaussian_elimination_partial_pivoting(A1,b1);

A2=[1 -1 2 -1;2 -2 3 -3;1 1 1 0;1 -1 4 3];
b2=[-8;-20;-2;4];
%x2=Gaussian_elimination_partial_pivoting(A2,b2)

A3=[pi -exp(1) sqrt(2) -sqrt(3);pi exp(1) -(exp(1)^2) 3/7;sqrt(5) -sqrt(6) 1 -sqrt(2);pi^3 exp(1)^2 -sqrt(7) 1/9];
b3=[sqrt(11);0;pi;sqrt(2)];
%x3=Gaussian_elimination_scaled_partial_pivoting(A3,b3)

x4=Gaussian_elimination_complete_partial_pivoting(A3,b3)

function x=Gaussian_elimination_partial_pivoting(A,b)
[p,n]=size(A);
if p~=n
    error('Not a square matrix; cannot compute soltuion');
end
if det(A)==0
    error('The determinant of the matrix is 0.');
end
A=[A b];
for i=1:1:n-1 % take ith row as reference
    % doing partial pivoting
    max=intmin;
    swap=0;
    for p=i:1:n
        if (max<abs(A(p,i)))
            max=abs(A(p,i));
            swap=p;
        end
    end
    if swap==0
        error("There is no unique solution");
    end
    if swap~=i
        temp=A(i,:);
        A(i,:)=A(swap,:);
        A(swap,:)=temp;
        fprintf("Step %i: Swapped row %i and row %i\n",i,i,swap);
    else
        fprintf("Step %i: No row swapped\n",i);
    end
    for j=i+1:1:n % perform operation on the jth row
        m=A(j,i)/A(i,i);
        for k=i:1:n+1
            A(j,k)=A(j,k)-m*A(i,k);
        end
    end
end
x=zeros(n,1);
x(n)=A(n,n+1)/A(n,n);
for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(A(i,n+1)-sum)/A(i,i);
end
A
end

    % if A(i,i)==0
    %     swap=0;
    %     for p=i+1:1:n
    %         if A(p,i)~=0
    %             swap=p;
    %             break;
    %         end
    %     end
    %     temp=A(i,:);
    %     A(i,:)=A(swap,:);
    %     A(swap,:)=temp;
    % end

function x=Gaussian_elimination_scaled_partial_pivoting(A,b)
[p,n]=size(A);
if p~=n
    error('Not a square matrix; cannot compute soltuion');
end
if det(A)==0
    error('The determinant of the matrix is 0.');
end
scale=max(abs(A'))'
A=[A b];
for i=1:1:n-1 % take ith row as reference
    % doing partial pivoting
    maximum=intmin;
    swap=0;
    for p=i:1:n
        if maximum<(abs(A(p,i))/scale(p))
            maximum=abs(A(p,i))/scale(p);
            swap=p;
        end
    end
    if swap==0
        error("There is no unique solution");
    end
    if swap~=i
        temp=A(i,:);
        A(i,:)=A(swap,:);
        A(swap,:)=temp;
        temp_scale=scale(i);
        scale(i)=scale(swap);
        scale(swap)=temp_scale;
        fprintf("Step %i: Swapped row %i and row %i\n",i,i,swap);
    else
        fprintf("Step %i: No row swapped\n",i);
    end
    for j=i+1:1:n % perform operation on the jth row
        m=A(j,i)/A(i,i);
        for k=i:1:n+1
            A(j,k)=A(j,k)-m*A(i,k);
        end
    end
end
x=zeros(n,1);
x(n)=A(n,n+1)/A(n,n);
for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(A(i,n+1)-sum)/A(i,i);
end
A
end

function x=Gaussian_elimination_complete_partial_pivoting(A,b)
[p,n]=size(A);
if p~=n
    error('Not a square matrix; cannot compute soltuion');
end
if det(A)==0
    error('The determinant of the matrix is 0.');
end
A=[A b];
swappings=zeros(1,2*n);
idx=1;
for i=1:1:n-1 % take ith row as reference
    % doing partial pivoting
    maximum=intmin;
    swap_x=0;
    swap_y=0;
    for g=i:1:n
        for h=i:1:n
            if maximum<abs(A(g,h))
                maximum=abs(A(g,h));
                swap_x=g;
                swap_y=h;
            end
        end
    end
    if swap_x~=i
        temp=A(i,:);
        A(i,:)=A(swap_x,:);
        A(swap_x,:)=temp;
        fprintf("Step %i: Swapped row %i and row %i\n",i,i,swap_x);
    else
        fprintf("Step %i: No row swapped\n",i);
    end
    if swap_y~=i
        temp=A(:,i);
        A(:,i)=A(:,swap_y);
        A(:,swap_y)=temp;
        swappings(idx)=i;
        swappings(idx+1)=swap_y;
        idx=idx+2;
        fprintf("Step %i: Swapped column %i and column %i\n",i,i,swap_y);
    else
        fprintf("Step %i: No column swapped\n",i);
    end
    for j=i+1:1:n % perform operation on the jth row
        m=A(j,i)/A(i,i);
        for k=i:1:n+1
            A(j,k)=A(j,k)-m*A(i,k);
        end
    end
end
x=zeros(n,1);
x(n)=A(n,n+1)/A(n,n);
for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(A(i,n+1)-sum)/A(i,i);
end
A
x
swappings
for q=idx-2:-2:1
    temp=x(swappings(q));
    x(swappings(q))=x(swappings(q+1));
    x(swappings(q+1))=temp;
end
end