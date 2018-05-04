function [bugs, bugsCellLabel] = CycIFBugSeg(FOVstack, cells)
%% bug segmentation
%creates a mask of bug location 
%but does not try to separate individual bugs
%CHANGE TO MAKE DEPENDENT ON CHANNEL NAME
B=FOVstack(:,:,3);

bugsBW=B>300;
bugs = bwlabel(bugsBW);
bugsCellLabel = cells.*uint16(bugsBW);

