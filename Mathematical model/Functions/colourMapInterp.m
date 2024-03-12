function [colourMap] = colourMapInterp(col1, col2, col3, n)
%GREYTOGREEN Summary of this function goes here
%   Detailed explanation goes here
% colourMap = [linspace(col1(1), col2(1),n)', linspace(col1(2), col2(2),n)', linspace(col1(3), col2(3),n)'];

colourMapFirstHalf = [linspace(col1(1), col2(1),n/2)', linspace(col1(2), col2(2),n/2)', linspace(col1(3), col2(3),n/2)'];
colourMapSecondHalf = [linspace(col2(1), col3(1),n/2)', linspace(col2(2), col3(2),n/2)', linspace(col2(3), col3(3),n/2)'];
colourMap = [colourMapFirstHalf; colourMapSecondHalf];
end

