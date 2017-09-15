function [PgradEb,Ebend]=BendingGrad(A,Kbend,ls)

%inputs,
   % A=n*2 array containing the n number of points, the initial MT guess
   % Kbend=relative weight of bending energy. Higher Kbend 
% for n=5 points (in m*2 array A), there are m=4 vectors between them
%A=[A(:,2),A(:,1)];
S=diff(A);

[m,~]=size(S);

D=zeros(1,m-1);
N=zeros(1,m-1);
% singleN=zeros(1,m);
% 
% for i=1:m
%     singleN(i)=norm(S(i,:));
% end
% L=sum(singleN);

Ebend = 0;
for i=1:(m-1)
D(i)=dot(S(i,:),S((i+1),:)); %dot products of adjacent vectors "i"
N(i)=norm(S(i,:))*norm(S((i+1),:)); %norms of adjacent vectors
Ebend = Ebend + Kbend/ls*(1-D(i)/N(i));
end

% Since there are m-1 terms in Ebend, there will be m-2 Ebend derivatives.
SdervEb=zeros(m-2,2);
for i=1:m-2
    for b=1:2
    SdervEb(i,b)=((S(i,b)/N(i)) + (S(i+2,b)/N(i+1)) - ((S(i+1,b)/(norm(S(i+1,:))^2)))*(D(i)/(N(i)) + D(i+1)/(N(i+1))));
    end
end


% Need to calculate the derivative (gradient) wrt the first and last S
% vectors separately
S1dervEb=zeros(1,2);
for b=1:2
S1dervEb(b)=(S(2,b)/N(1)) - S(1,b)/norm(S(1,:))^2*D(1)/N(1);
end

SmdervEb=zeros(1,2);
for b=1:2
SmdervEb(b)=(S(m-1,b)/N(m-1)) - S(m,b)/norm(S(m,:))^2*D(m-1)/N(m-1);
end

SgradEb=zeros(m-1,2);
SgradEb=[S1dervEb; SdervEb; SmdervEb]; %gradient of bending energy wrt vectors S

%Need gradient wrt Cartesian coordinates i.e. the m number of points (x,y)
PdervEb=zeros(m-2,2);
for i=1:m-1
    for b=1:2
    PdervEb(i,b)=SgradEb(i,b)-SgradEb(i+1,b);
    end
end

% Need to calculate the derivative wrt the first and last points separately
P1dervEb=-S1dervEb;
PmdervEb=SmdervEb;

%Gradient wrt Cartesian coordinates (discrete pixels)
PgradEb=-(Kbend/ls)*[P1dervEb; PdervEb; PmdervEb];

% I=grayscale/RGB image
% Quiver generates arrows at each of the points in array A guiding snake
% toward microtubules

% if (nargin>2)
% imshow(I), hold on; quiver(A(:,2),A(:,1),PgradEb(:,2),PgradEb(:,1),'m')
% end