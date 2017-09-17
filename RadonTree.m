function [angle2compare] = RadonTree(img, pstart, accIntensity,  position, angle2compare, opt)

global path_vertexs;
global max_pos;
global max_compare
%%
% opt = struct();
% opt.Lradon = 70; % length to use for finite radon transform
% opt.Ltrace = 5; % separation between points for the tracing
% % cutoff for radon transform magnitude to call filament end
% % relative to starting point magnitude
% % cutoff in theta change to call filament end
% % in degrees
% opt.relthcutoff = 60;
% 
% % max points to trace in each direction
% opt.maxpts = 100;
% 
% % display tracing during calculation
% opt.dodisplay = 1;




img3d = double(img);

L = opt.Lradon; % 100

[validated_peak_angle] = RadonTransform(img3d,pstart,L, angle2compare, opt);

if (validated_peak_angle==0)
    [compare] =  pathEvaluate(accIntensity, position)
    %end values and sum and go to new function => evaluator
    %compare values and find out max value - then plot the positions and
    %display brightness/energy values
     if compare > max_compare
         max_compare = compare
         max_pos = position;
     end
    
else
    for i = 1:size(validated_peak_angle,2)
%         size(validated_peak_angle)
        accIntensity = accIntensity + validated_peak_angle(2,i);
        %update pnext using pValue
        
    
        theta = validated_peak_angle(1,i)*pi/180;
        dvec = [cos(theta),sin(theta)];
        
        
        if (opt.dodisplay)
            lams = 0:1:opt.Ltrace;
            linepts = [dvec(1)*lams'+pstart(1),dvec(2)*lams'+pstart(2)];
            hold all
%             plot(pstart(1),pstart(2),'R*')
%             pause (1)
            plot(pstart(1),pstart(2),'k*')
%             figure('visible','off')
%             pause (1)
%             figure('visible','on')
            plot(pstart(1),pstart(2),'y.',linepts(:,1),linepts(:,2),'g')
            plot(linepts(end,1),linepts(end,2),'R*')
            hold off
            drawnow
        end
        
        pnext = pstart + dvec*opt.Ltrace;
        d1 = path_vertexs - repmat(pstart,size(path_vertexs,1),1);
        d2 = path_vertexs - repmat(pnext,size(path_vertexs,1),1);
        d1 = sum(d1.^2,2);
        d2 = sum(d2.^2,2);
        r1=min(d1); 
        r2=min(d2); 
        if r1 < opt.path_residue && r2<opt.path_residue
            return;
        end
        position = [position; pnext];
        path_vertexs = [path_vertexs; pnext];
        RadonTree(img3d, pnext, accIntensity, position, validated_peak_angle(1,i), opt);
    end

end




