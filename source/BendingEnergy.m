%% Bending Energy
function Eb = BendingEnergy(A,Kbend,ls)
% A is a 2-column matrix containing all the coordinates corresponding to
% a microtubule.
% ls = ground state length of segment

% Example pts = [50 50; 51 47; 55 40; 60 41]

S=diff(A);
[m,~]=size(S);

D=zeros(1,m-1);
N=zeros(1,m-1);
Ebend=zeros(1,m-1);

% cos(x) = U*V/(u*v) where u = norm(U)
% Assume zero intrinsic curvature

for i=1:(m-1)
D(i)=dot(S(i,:),S((i+1),:)); %dot products of adjacent vectors "i"
N(i)=norm(S(i,:))*norm(S((i+1),:)); %norms of adjacent vectors
Ebend(i)=1 - D(i)/N(i); %removed prefactor Kbend*m/L for brevity
end

% singleN=zeros(1,m);
% for i=1:m
%     singleN(i)=norm(S(i,:));
% end
% L=sum(singleN); % Calculates the displacement between MT points

%Kbend is a constant with dimensions of energy/length(in pixels).
Eb=(Kbend/ls)*sum(Ebend);