%%
% generate a snake
pts = [0 0; 1 1; 2 0; 3 1];
pts = pts + rand(size(pts))*1;

plot(pts(:,1),pts(:,2),'.-')

%%
pts = [0 0; 1 1; 2 2];
plot(pts(:,1),pts(:,2),'.-')
%%
Kbend = 1; ls = 1;
%Eb0 = BendingEnergy(pts,Kbend,ls);

[PgradEb0, Eb0] =BendingGrad(pts,Kbend,ls);
%
tiny = 1e-6;
for pc= 1:size(pts,1)
    for xc = 1:2
        pts(pc,xc) = pts(pc,xc) + tiny;
        [PgradEb, Eb] =BendingGrad(pts,Kbend,ls);
        pts(pc,xc) = pts(pc,xc) - tiny;
        
        [pc xc (Eb-Eb0)/tiny PgradEb0(pc,xc)]
    end
end

% ----------------------------


% -----------------
%% function for getting derivatives using chain object

chain = struct();
chain.ls = sqrt(2); % preferred segment length
chain.lp = 1; % bending persistence length
chain.lstretch = 1; % stretching energy prefactor
% don't fix end orientations
chain.fixtan0 = 0; 
chain.fixtanf = 0;
chain.tan0 = zeros(2,1);
chain.tanf = zeros(2,1);
% WARNING: currently fixes end positions
chain.pos0 = pts(1,:)';
chain.posf = pts(end,:)';

coords = reshape(pts(2:end-1,:)',(size(pts,1)-2)*2,1);
chain.ncrd = length(coords);
chain.nseg = size(pts,1)-1;
chain.nbead = size(pts,1)-2;
[energy,grad] = energyWLC2d(chain,coords)

%% check derivatives

[energy0,grad0] = energyWLC2d(chain,coords)
%
tiny = 1e-4;
for pc= 1:chain.ncrd
        coords(pc) = coords(pc)+tiny;
        [energy] = energyWLC2d(chain,coords);
        coords(pc) = coords(pc)-tiny;       
        
        [pc (energy-energy0)/tiny grad0(pc,xc)]
end

%%
grad = reshape(grad0,chain.nbead,2)';

plot(pts(:,1),pts(:,2),'.-')
hold all
quiver(pts(2:end-1,1),pts(2:end-1,2),grad(:,1),grad(:,2))
hold off