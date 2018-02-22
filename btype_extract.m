function [IDout,btype] = btype_extract(ID)
% extract the btype id (0,1) and matrix case ID from the two digit number
% [ID,btype]

btype = mod(ID,10);
IDout = (ID - btype)/10;