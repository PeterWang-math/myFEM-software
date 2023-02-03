function [ufem, outindex] = utriFEM2(mesh,basis,uh,x,y,s,t,femtype,upperbd)
% generate FEM solution function with s+t order derevatives, i.e. 
% y = d^{s+t}(uh)/dx^{s}dy^{t} 
% x & y can be Mx1 column vector
% outindex returns the points' index that are not in \Omega
% upperbd: maximum amount in each batch

ufem = 0;
M = length(x); % number of points

trnlb = basis.trial.Nlb;
trtb = basis.trial.Tb;
trf = basis.trial.reffun;

% get element number
if M < upperbd
    Nindex = pt2elem([x,y],mesh.P,mesh.T);
else
    subnum = floor(M/upperbd);
    for i = 1:subnum
        Nindex(upperbd*(i-1)+1:upperbd*i) = pt2elem([x(upperbd*(i-1)+1:upperbd*i),...
                                            y(upperbd*(i-1)+1:upperbd*i)],mesh.P,mesh.T);
    end
    if subnum*upperbd<M
        Nindex(subnum*upperbd+1:M) = pt2elem([x(subnum*upperbd+1:end),...
                                       y(subnum*upperbd+1:end)],mesh.P,mesh.T);
    end
end
[~, outindex] = find(Nindex == -1); % points that are not in \Omega 
Nindex(outindex) = 1; % insure the legal input 

% generate block diagonal matrix
v21 = mesh.P(:,mesh.T(2,Nindex)) - mesh.P(:,mesh.T(1,Nindex));
v31 = mesh.P(:,mesh.T(3,Nindex)) - mesh.P(:,mesh.T(1,Nindex));
detJ = v21(1,:).*v31(2,:) - v31(1,:).*v21(2,:);
J = [reshape(v21,[],1),reshape(v31,[],1)];
v01 = J(1:2:end-1,1)';
v02 = J(2:2:end,2)';
v0 = reshape([v01;v02],[],1);
v1 = reshape([J(1:2:end-1,2)';zeros(1,M)],[],1);
vn1 = reshape([J(2:2:end,1)';zeros(1,M)],[],1);
%J = diag(v0,0) + diag(v1(1:end-1),1) + diag(vn1(1:end-1),-1);
J = sparse([2:2*M,1:2*M,1:2*M-1],[1:2*M-1,1:2*M,2:2*M],[vn1(1:end-1)',v0',v1(1:end-1)']);

% compute phat
p = reshape([x,y]'-mesh.P(:,mesh.T(1,Nindex)),[],1);
phat = reshape(J\p,[],M)';

s = num2str(s);  t = num2str(t); 
switch femtype.type
    case 'pgm'
        switch strcat(s,t)
            case '00'
                U = zeros(M,trnlb);
                for k = 1:trnlb
                    U(:,k) = trf{k,1,1}(phat(:,1),phat(:,2));
                end
                ufem = sum( U.*uh(trtb(:,Nindex))', 2);  
            case '10'
                J22 = v31(2,:)./detJ;
                nJ21 = -v21(2,:)./detJ;
                Uxhat = zeros(M,trnlb);
                Uyhat = Uxhat;
                for k = 1:trnlb
                    Uxhat(:,k) = trf{k,2,2}(phat(:,1),phat(:,2));
                    Uyhat(:,k) = trf{k,1,2}(phat(:,1),phat(:,2));
                end
                ufem = sum( Uxhat.*( repmat(J22',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uyhat.*( repmat(nJ21',1,trnlb).*uh(trtb(:,Nindex))' ), 2);

            case '01'
                nJ12 = -v31(1,:)./detJ;
                J11 = v21(1,:)./detJ;
                Uxhat = zeros(M,trnlb);
                Uyhat = Uxhat;
                for k = 1:trnlb
                    Uxhat(:,k) = trf{k,2,2}(phat(:,1),phat(:,2));
                    Uyhat(:,k) = trf{k,1,2}(phat(:,1),phat(:,2));
                end
                ufem = sum( Uxhat.*( repmat(nJ12',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uyhat.*( repmat(J11',1,trnlb).*uh(trtb(:,Nindex))' ), 2) ;
                
            case '20'
                J22 = v31(2,:)./detJ;
                nJ21 = -v21(2,:)./detJ;
                cn1 = J22.^2;  cn2 = 2*J22.*nJ21;  cn3 = nJ21.^2;
                Uxxhat = zeros(M,trnlb);
                Uxyhat = Uxxhat;
                Uyyhat = Uxxhat;
                for k = 1:trnlb
                    Uxxhat(:,k) = trf{k,3,3}(phat(:,1),phat(:,2));
                    Uxyhat(:,k) = trf{k,2,3}(phat(:,1),phat(:,2));
                    Uyyhat(:,k) = trf{k,1,3}(phat(:,1),phat(:,2));
                end
                ufem = sum( Uxxhat.*( repmat(cn1',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uxyhat.*( repmat(cn2',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uyyhat.*( repmat(cn3',1,trnlb).*uh(trtb(:,Nindex))' ), 2);

            case '11'
                J22 = v31(2,:)./detJ;
                nJ21 = -v21(2,:)./detJ;
                nJ12 = -v31(1,:)./detJ;
                J11 = v21(1,:)./detJ;
                cn1 = nJ12.*J22;  cn2 = nJ12.*nJ21 + J11.*J22;  cn3 = J11.*nJ21;
                Uxxhat = zeros(M,trnlb);
                Uxyhat = Uxxhat;
                Uyyhat = Uxxhat;
                for k = 1:trnlb
                    Uxxhat(:,k) = trf{k,3,3}(phat(:,1),phat(:,2));
                    Uxyhat(:,k) = trf{k,2,3}(phat(:,1),phat(:,2));
                    Uyyhat(:,k) = trf{k,1,3}(phat(:,1),phat(:,2));
                end
                ufem = sum( Uxxhat.*( repmat(cn1',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uxyhat.*( repmat(cn2',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uyyhat.*( repmat(cn3',1,trnlb).*uh(trtb(:,Nindex))' ), 2);

            case '02'
                nJ12 = -v31(1,:)./detJ;
                J11 = v21(1,:)./detJ;
                cn1 = nJ12.^2;  cn2 = 2*nJ12.*J11;  cn3 = J11.^2;
                Uxxhat = zeros(M,trnlb);
                Uxyhat = Uxxhat;
                Uyyhat = Uxxhat;
                for k = 1:trnlb
                    Uxxhat(:,k) = trf{k,3,3}(phat(:,1),phat(:,2));
                    Uxyhat(:,k) = trf{k,2,3}(phat(:,1),phat(:,2));
                    Uyyhat(:,k) = trf{k,1,3}(phat(:,1),phat(:,2));
                end
                ufem = sum( Uxxhat.*( repmat(cn1',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uxyhat.*( repmat(cn2',1,trnlb).*uh(trtb(:,Nindex))' ) + ...
                            Uyyhat.*( repmat(cn3',1,trnlb).*uh(trtb(:,Nindex))' ), 2);
        end
    case 'ipdg'
        
        
end
ufem(outindex) = nan;
end


function elemnum = pt2elem(pt,P,T)
% subroutine: given a mesh and a point (x,y),find which element the point belongs to 
% elemnum (1xM) returns the result, if \Omega doesn't contain pt, elemnum = -1;
% pt is a 1x2 point or Mx2 vector
% P: mesh point
% T: mesh element

M = size(pt,1); % number of points
N = size(T,2); % number of elements
elemnum = zeros(1,M);
%% vectorization programming, takes up a lot of memory
pt = reshape( repmat(pt',N,1), [], N*M );

s11 = repmat(P(:,T(1,:)),1,M) - pt;
s12 = repmat(P(:,T(2,:)),1,M) - pt;
S1 = s11(1,:).*s12(2,:) - s12(1,:).*s11(2,:);
clear s11 s12
s22 = repmat(P(:,T(2,:)),1,M) - pt;
s23 = repmat(P(:,T(3,:)),1,M) - pt;
S2 = s22(1,:).*s23(2,:) - s23(1,:).*s22(2,:);
clear s22 s23
s33 = repmat(P(:,T(3,:)),1,M) - pt;
s31 = repmat(P(:,T(1,:)),1,M) - pt;
S3 = s33(1,:).*s31(2,:) - s31(1,:).*s33(2,:);
clear s33 s31 pt 

index = S1>=0 & S2>=0 & S3>=0;
index = reshape(index',[],M);
clear S1 S2 S3

outnum = find(sum(index,1)==0);
if isempty(outnum) % all points are in \Omega
    [~,I] = max(index); % find first non-zero index in each column as elemnum
    elemnum = I;
else
    leftindex = 1:M; leftindex(outnum) = [];
    elemnum(outnum) = -1; % points that are not in \Omega
    [~,I] = max(index(:,leftindex)); % find first non-zero index in each column as elemnum
    elemnum(leftindex) = I;
end

end


