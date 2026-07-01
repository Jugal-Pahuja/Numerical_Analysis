y1=system_elements(2,3,-1,2)
plot_system(2,3,-1,2)
y2=system_elements(10,1,1,1)
plot_system(10,1,1,1)

function y=UFL(beta,L)
y=beta^(L-1);
end

function y=OFL(beta,t,U)
y=(1-(1/(beta^t)))*(beta^U);
end

function y=no_of_elements(beta,t,L,U)
y=2*(beta-1)*(beta^(t-1))*(U-L+1)+1;
end

function y=system_elements(beta,t,L,U)
% the total number of elements in the system
n=2*(beta-1)*(beta^(t-1))*(U-L+1)+1;
% vector which stores all the elements
y=zeros(n,1); 
% forming the matrix which makes the mantissa showing d2 to dn
matrix=zeros(beta^(t-1),t-1);
for k=1:t-1 % for making each coloumn
    for times=0:beta^(t-k-1)-1
        for i=0:beta-1
            for j=1:beta^(k-1)
                matrix(times*(beta^k)+i*(beta^(k-1))+j,k)=i;
            end
        end
    end
end
index=2;
for e=L:U
    for i=1:beta^(t-1)
        for j=1:beta-1
            s=0;
            for k=1:t-1
                s=s+matrix(i,k)/(beta^(k+1));
            end
            s=s+j/beta;
            s=s*(beta^e);
            y(index,1)=s;
            y(index+1,1)=-1*s;
            index=index+2;
        end
    end
end
y=sort(y);
end

function []=plot_system(beta,t,L,U)
figure;
y=system_elements(beta,t,L,U);
plot(y,0,'b.');
xlabel('System Elements');
title('Plot of all elements of F(beta,t,L,U)');
end

