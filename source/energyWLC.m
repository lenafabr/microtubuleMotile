function [energy,grad] = energyWLC(chain,coords)
% get the energy and gradient for a wormlike chain
% copied from 
%/storageserver/backups/SpakowitzBackups/phoenix/Spakowitz20150907/Spakowitz/proj/rnap/optimize
% modified to use fixed end positions and tangents at both ends
% starting position / tangent: chain.pos0, chain.tan0
% ending position / tangent: chain.pof, chain.tanf


if nargin == 1
    coords = chain.coords;
end


LP = chain.lp/chain.ls;
LST = chain.lstretch/chain.ls/2;

Estretch = 0; Ebend = 0;
Gstretch = zeros(chain.ncrd,1); Gbend = zeros(chain.ncrd,1);

prevseg = chain.tan0; nprev = 1;

% ---------- (-) end of chain -----------------
% initial position and tangent fixed to chain.pos0, chain.tan0
if (chain.nseg > 1)
    % get the stretch energy for the first segment
    seg = coords(1:3)-chain.pos0; nseg = norm(seg);
    diff = nseg-chain.ls;
    Estretch = Estretch + LST*diff^2;
    Gstretch(1:3) = Gstretch(1:3) + 2*LST*diff*seg/nseg;
    
    if (chain.fixtan0)
        % and the bend energy at the fixed bead
        cang = (seg'*prevseg)/nseg/nprev;
        Ebend = Ebend + LP*(1-cang);
        % derivatives wrt elements of seg
        Gbend(1:3) = Gbend(1:3)-LP*(prevseg/nprev - cang*seg/nseg)/nseg;
    end
    
    prevseg = seg; nprev = nseg;
end
% -------------------------------------------------

% --------------------
% get the energy and gradient for the internal beads
% ---------------------

for bc = 1:chain.nbead
    if (bc==chain.nbead)
        seg = chain.posf - coords(end-2:end);
    else
        seg = coords(3*bc+1:3*(bc+1))-coords(3*(bc-1)+1:3*bc);
    end
    
    % calculate stretch energy from bc->bc+1 segment
    nseg = norm(seg);
    diff = nseg-chain.ls;
    Estretch = Estretch+LST*diff^2;
    % stretch gradient wrt seg
    addgrad = 2*LST*diff*seg/nseg;
    if (bc < chain.nbead)
        Gstretch(3*bc+1:3*(bc+1)) = Gstretch(3*bc+1:3*(bc+1)) + addgrad;
    end
    Gstretch(3*(bc-1)+1:3*bc) = Gstretch(3*(bc-1)+1:3*bc) - addgrad;
    
    % calculate bend energy at bc-1 -> bc -> bc+1 junction
    % cosine of angle between segments
    cang = (seg'*prevseg)/nseg/nprev;
    Ebend = Ebend + LP*(1-cang);
    % bend gradient wrt seg
    gradS = -LP*(prevseg/nprev - cang*seg/nseg)/nseg;
    % bend gradient wrt prevseg
    gradP = -LP*(seg/nseg - cang*prevseg/nprev)/nprev;
    
    if (bc<chain.nbead)
        Gstretch(3*bc+1:3*(bc+1)) = Gstretch(3*bc+1:3*(bc+1)) + gradS;
    end
    Gstretch(3*(bc-1)+1:3*bc) = Gstretch(3*(bc-1)+1:3*bc) - gradS + gradP;
    if (bc>1) 
        Gstretch(3*(bc-2)+1:3*(bc-1)) = Gstretch(3*(bc-2)+1:3*(bc-1))-gradP;
    end
    
    prevseg = seg; nprev = nseg;   
end

% ---------- (+) end of chain -----------------
% final position and tangent fixed to chain.posf, chain.tanf
if (chain.fixtanf && chain.nseg > 1)    
    % and the bend energy at the fixed bead
    ntanf = norm(chain.tanf);
    cang = (chain.tanf'*prevseg)/ntanf/nprev;
    Ebend = Ebend + LP*(1-cang);
    % derivatives wrt elements of seg
    Gbend(end-2:end) = Gbend(end-2:end)+LP*(chain.tanf/ntanf - cang*prevseg/nprev)/nprev;    
end
% -------------------------------------------------


energy = Estretch+Ebend;
grad = Gstretch+Gbend;

%[energy norm(grad)]

%[energy,norm(grad)]
end