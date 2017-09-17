%% Bending Energy
function Eb = BendingEnergy(A)
% A is a 2-column matrix containing all the coordinates corresponding to
% a microtubule.

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
Ebend(i)=D(i)/N(i)+1; %removed prefactor Kbend*(m-1)/L for brevity
end

% Sets one of the terminal coordinates as the origin.
L=norm(S); % Calculates the displacement between MT end points


%Kbend is a constant with dimensions of energy/length(in pixels).
Eb = ((m-1)/L)*sum(Ebend+1);