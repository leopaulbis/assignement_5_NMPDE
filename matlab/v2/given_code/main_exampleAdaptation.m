clear all; close all;
addpath('src_adapt/')

disp('The objective of this source code is to illustrate the use of refineListElements')

doPlot = 1;

figure(1);
[mesh]=generateMeshSquare(doPlot);

numRefs = 1;

figure(2);
for iref=1:numRefs
    
    % For illustration purposes:
    %   we select to refine elements close to a circle
    cm   = (mesh.X(mesh.T(:,1),:)+mesh.X(mesh.T(:,2),:)+mesh.X(mesh.T(:,3),:))/3.0;
    R    = 0.25;
    r    = sqrt(sum((cm-0.5).^2,2));
    tolr = 0.5/(2^iref);
    listElems = find(abs(r-R)<tolr);
    disp(listElems);
    [mesh]=refineListElements(mesh,listElems,doPlot);  
    
    pause(1)
end

disp('How to use mesh data structure:')
disp('  mesh.X  -> coordiantes')
disp('  mesh.T  -> connectivities')
disp('  mesh.boundaryNodes  -> boundary nodes id')
