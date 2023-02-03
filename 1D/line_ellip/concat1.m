function [row,col,s] = concat1(row1,col1,s1,row2,col2,s2,trnlb,tsnlb,ordtr,ordts)
% concat two array s1,s2 to s
if isempty(row1)==1
    row = row2;
    col = col2;
    s = s2;
    return
end
l1 = length(row1);
l2 = length(row2);

% get the index of overlapping block between two elements
% block1 = reshape(l1-tsnlb*ordtr+1:l1,tsnlb,ordtr);
% iamat = block1(tsnlb-ordts+1:end,1:end);
% x1mat = block1(1:tsnlb-ordts,1:ordtr);
% ia = reshape(iamat,1,ordts*ordtr);
% 
% block2 = reshape(1:ordtr*tsnlb,tsnlb,ordtr);
% ibmat = block2(1:ordts,1:ordtr);
% x2mat = block2(ordts+1:tsnlb,1:ordtr);
% ib = reshape(ibmat,1,ordts*ordtr);

X1 = row1+col1*1i;
X2 = row2+col2*1i;
[~,ia,ib] = intersect(X1,X2);

l3 = length(ia); 
s = zeros(1,l1+l2-l3);

% get the remaining index of s1 and s2
x1 = setdiff(1:l1,ia);
x2 = setdiff(1:l2,ib);
% x1 = [1:l1-ordtr*tsnlb,reshape(x1mat,1,ordtr*(tsnlb-ordts))];
% x2 = [reshape(x2mat,1,(tsnlb-ordts)*ordtr),tsnlb*ordtr+1:l2];

row(1:l1-l3) = row1(x1);
row(l1-l3+1:l1) = row1(ia);
row(l1+1:l1+l2-l3) = row2(x2);

col(1:l1-l3) = col1(x1);
col(l1-l3+1:l1) = col1(ia);
col(l1+1:l1+l2-l3) = col2(x2);

s(1:l1-l3) = s1(x1);
s(l1-l3+1:l1) = s1(ia) + s2(ib);
s(l1+1:l1+l2-l3) = s2(x2);

end

