close all
clc
format loose
n=16;            % Order of the graph

B= [[0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]];

A= reshape(B,[n,n]);
if A==A'
      disp('Matrix is Symmetric')
else
      disp('Not Symmetric')
end

l=size(A,1);
e=[];
s=[];
t=[];
k=[];
for i=1:l
    for j=i+1:l
        if A(i,j)==1
            e=[e;[i,j]];
        end
    end
end

edges=size(e,1);
j=ones(edges,1);

k=sum(A);
k1=k';

p=[];
for i=1:edges
    pp=[e(i,:) k(e(i,1)) k(e(i,2))];
    p=[p;pp];
end

A1=p(:,3);
B1=p(:,4);
A2= A1.*B1;
D2= A1.^2 + B1.^2;
E2= sqrt(D2); 
C2= (A1-1).*(B1-1);
B2=A1+B1;

A3=B2-2*j;
A4=sqrt(A3./A2);
ABC1=sum(A4);                               % ABC index
% 
A5=(A2./A3).^3;
AZI=sum(A5);                                % AZI index
% 
A6=(2.*sqrt(A2))./B2;
GA=sum(A6);                                  % GA index
% 
A7=A2.^(-1);
R1=sum(A7);                                 % R(-1) index
% 

A9=1/sqrt(A2);                                % Randic index
R=sum(A9);

A10=sum(A2);                                % R_{1} index

A11=(A2).^2;                                % R_{2} index
R2=sum(A11);

A12=(A2).^(-2);                             % R_{-2} index
R21=sum(A12);

A13=(A2).^(-0.2661);                        %R_{-0.2661} index
R22=sum(A13);

%
A8=B2;
SCI1=sum(A8);                               % SCI(1) index
A14=B2.^(2);
SCI2=sum(A14);                               % SCI(2) index
A15=B2.^(-1);
SCI12=sum(A15);                               % SCI(-1) index
A16=B2.^(-2);
SCI13=sum(A16);                               % SCI(-2) index
A17=1./(B2.^(1/2));
SCI14=sum(A17);                               % SCI(-1/2) index
A18=B2.^(0.5601);
SCI15=sum(A18);                               % SCI(-0.5601) index
%
A19=A2.^(1/2);                                  %RR(G) index
RR=sum(A19);
%
A20=C2.^(1/2);                                  %RRR(G) index
RRR=sum(A20);
%
A21=B2;                                         %M_{1}(G) index
M1=sum(A21);
A22=A2;                                         %M_{2}(G) index
M2=sum(A22);
%
A23=B2;                                         %Pi_1(G) Index
Pi1= prod(A23);
A24=A2;                                         %Pi_2(G) Index
Pi2= prod(A24);
%
SO= sum(E2);                                    %Somber index

disp('Degree-based indices:')
fprintf('The Randic Index is %4.4f\n',R');
fprintf('The R_{1} Index is %4.4f\n',A10');
fprintf('The R_{-1} Index is %4.4f\n',R1');
fprintf('The R_{2} Index is %4.4f\n',R2');
fprintf('The R_{-2} Index is %4.4f\n',R21');
fprintf('The R_{-0.2661} Index is %4.4f\n',R22');
fprintf('The SCI_{1} Index is %4.4f\n',SCI1');
fprintf('The SCI_{-1} Index is %4.4f\n',SCI12');
fprintf('The SCI_{2} Index is %4.4f\n',SCI2');
fprintf('The SCI_{-2} Index is %4.4f\n',SCI13');
fprintf('The SCI_{-1/2} Index is %4.4f\n',SCI14');
fprintf('The SCI_{-0.5601} Index is %4.4f\n',SCI15');
fprintf('The ABC Index is %4.4f\n',ABC1');
fprintf('The AZI Index is %4.4f\n',AZI');
fprintf('The GA Index is %4.4f\n',GA');
fprintf('The RR(G) Index is %4.4f\n',RR');
fprintf('The RRR(G) Index is %4.4f\n',RRR');
fprintf('The M_{1}(G) Index is %4.4f\n',M1');
fprintf('The M_{2}(G) Index is %4.4f\n',M2');
fprintf('The Pi_{1}(G) Index is %4.4f\n',Pi1');
fprintf('The Pi_{2}(G) Index is %4.4f\n',Pi2');
fprintf('The SO(G) Index is %4.4f\n',SO');