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

if (chain.fixtan0)
    prevseg = chain.tan0; 
end
nprev = 1;

% ---------- (-) end of chain -----------------
% initial position and tangent fixed to chain.pos0, chain.tan0
if (chain.nseg > 1)
    % get the stretch energy for the first segment
    if (chain.fixpos0)
        seg = coords(1:2)-chain.pos0; nseg = norm(seg);
        diff = nseg-chain.ls;
        Estretch = Estretch + LST*diff^2;
        Gstretch(1:2) = Gstretch(1:2) + 2*LST*diff*seg/nseg;
        prevseg = seg; nprev = nseg;
    end
    
    if (chain.fixtan0)
        if (~chain.fixpos0)
            error('not set up for fixed gradients but not positions')
        end
        % and the bend energy at the fixed bead
        cang = (seg'*prevseg)/nseg/nprev;
        Ebend = Ebend + LP*(1-cang);
        % derivatives wrt elements of seg
        Gbend(1:2) = Gbend(1:2)-LP*(prevseg/nprev - cang*seg/nseg)/nseg;
    end
    
end
% -------------------------------------------------

% --------------------
% get the energy and gradient for the internal beads
% ---------------------

for bc = 1:chain.nbead
    if (bc==chain.nbead)
        if (chain.fixposf)
            seg = chain.posf - coords(end-1:end);
        else
            break
        end
    else
        seg = coords(2*bc+1:2*(bc+1))-coords(2*(bc-1)+1:2*bc);
    end
    
    % calculate stretch energy from bc->bc+1 segment
    nseg = norm(seg);
    diff = nseg-chain.ls;
    Estretch = Estretch+LST*diff^2;
    % stretch gradient wrt seg
    addgrad = 2*LST*diff*seg/nseg;
    if (bc < chain.nbead)
        Gstretch(2*bc+1:2*(bc+1)) = Gstretch(2*bc+1:2*(bc+1)) + addgrad;
    end
    Gstretch(2*(bc-1)+1:2*bc) = Gstretch(2*(bc-1)+1:2*bc) - addgrad;
    
    % calculate bend energy at bc-1 -> bc -> bc+1 junction
    % cosine of angle between segments
    if (~chain.fixpos0 & bc==1)
        prevseg = seg; nprev = nseg;
       continue
    end
    cang = (seg'*prevseg)/nseg/nprev;
    Ebend = Ebend + LP*(1-cang);
    % bend gradient wrt seg
    gradS = -LP*(prevseg/nprev - cang*seg/nseg)/nseg;
    % bend gradient wrt prevseg
    gradP = -LP*(seg/nseg - cang*prevseg/nprev)/nprev;
    
    if (bc<chain.nbead)
        Gstretch(2*bc+1:2*(bc+1)) = Gstretch(2*bc+1:2*(bc+1)) + gradS;
    end
    Gstretch(2*(bc-1)+1:2*bc) = Gstretch(2*(bc-1)+1:2*bc) - gradS + gradP;
    if (bc>1) 
        Gstretch(2*(bc-2)+1:2*(bc-1)) = Gstretch(2*(bc-2)+1:2*(bc-1))-gradP;
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
    Gbend(end-1:end) = Gbend(end-1:end)+LP*(chain.tanf/ntanf - cang*prevseg/nprev)/nprev;    
end
% -------------------------------------------------


energy = Estretch+Ebend;
grad = Gstretch+Gbend;

%[energy norm(grad)]

%[energy,norm(grad)]
end