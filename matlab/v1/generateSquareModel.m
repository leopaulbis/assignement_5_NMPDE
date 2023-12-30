function [g]=generateSquareModel()

model = createpde;
R1 = [3,4,0,1,1,0,1,1,0,0]';
g = decsg(R1);
% geometryFromEdges(model,g);
% pdegplot(model,EdgeLabels="on")
% axis equal