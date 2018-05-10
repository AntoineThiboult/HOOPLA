function [MergedStruct]=mergeStruct(Struct1, Struct2)
%
% [MergedStruct]=mergeStruct(Struct1, Struct2)
% 
% Merge the two structure Struct1 and Struct2. If the two structures have
% fieldnames in common, the values of Struct2 replace the ones of Struct1
%

M = [fieldnames(Struct1)' fieldnames(Struct2)'; struct2cell(Struct1)' struct2cell(Struct2)'];
[~, rows] = unique(M(1,:), 'last');
M=M(:, rows);
MergedStruct=struct(M{:});
