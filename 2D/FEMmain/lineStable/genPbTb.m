function PT = genPbTb(mesh,type,order,femtype)
% generate FEM nodes infomation in 2d
%
% The following input infomation is needed:
% mesh: a struct consisting of mesh.P, mesh.T, mesh.N;
% type: type of FEM basis, the following type is considered:
%       'l' for Lagrange basis funciton,
%       'h' for Hermite basis function,
%       'cr' for Crouzeix-Raviart element,
%       'a' for Argyris element,
%       'b' for Bell element,
%       'hct' for Hsich-Clough-Tocher element,
%       'rt' for Raviart-Thomas element
% order: a positive integer representing the order of basis function, less than 8
% femtype: 'pgm': Petrov-Galerkin method; 
%          'ipdg': internal penalty discontinous Galerkin method;
%
% The output information includes:
% PT.Pb & PT.Tb: FEM nodes information

elemtype = size(mesh.T,1);
nodeNum = size(mesh.P,2);

%% generate Pb and Tb matrix, Tb should consist with local basis function
if elemtype == 3  % triangular mesh
    if contains(['l','h','cr','a','b','hct','rt'],type) == 0
        warning(['The type for Pb & Tb must be either l, h, cr, a, b, hct, or rt ',...
            'already reset to l']);
        type = 'l';
    end
    switch type
        case 'l'  % Lagrange type FEM element
            if floor(order)~=order || order<1 || order>7
                warning(['The order for Pb & Tb must be a positive integer less than 8,',...
                    'already reset to 1']);
                order = 1;
            end
            switch order
                case 1 % Lagrange linear basis function
                    switch femtype.type
                        case 'pgm'
                            PT.Pb = mesh.P;
                            PT.Tb = mesh.T;
                        case 'ipdg'
                            
                    end 
                case 2 % Lagrange quadratic basis function
                    switch femtype.type
                        case 'pgm'
                            Pb = mesh.P; Tb = [mesh.T ; mesh.T]; 
                            Tb(4,:) = Tb(2,:); Tb(6,:) = Tb(3,:); 
                            % generate all midnode
                            midnode14 = 1/2*(Pb(:,Tb(1,:))+Pb(:,Tb(4,:)));
                            midnode16 = 1/2*(Pb(:,Tb(1,:))+Pb(:,Tb(6,:)));
                            midnode46 = 1/2*(Pb(:,Tb(4,:))+Pb(:,Tb(6,:)));
                            midnode14 = midnode14(1,:)+1i*midnode14(2,:);
                            midnode16 = midnode16(1,:)+1i*midnode16(2,:);
                            midnode46 = midnode46(1,:)+1i*midnode46(2,:);
                            midnode = unique([midnode14, midnode16, midnode46],'stable');
                            Pb = [Pb, [real(midnode); imag(midnode)]];
                            % add midnode to Tb
                            [~,index1]=ismember(midnode14,midnode);
                            [~,index2]=ismember(midnode16,midnode);
                            [~,index3]=ismember(midnode46,midnode);
                            Tb(2,:) = index1 + nodeNum;
                            Tb(3,:) = index2 + nodeNum;
                            Tb(5,:) = index3 + nodeNum;
                            PT.Pb = Pb;
                            PT.Tb = Tb;
                        case 'ipdg'
                            
                    end
                case 3 % Lagrange cubic basis function
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 4 % Lagrange fourth-order basis function
                    switch femtype.type
                        case 'pgm'
                            
                            
                        case 'ipdg'
                            
                    end
                case 5 % Lagrange fifth-order basis function
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 6 % Lagrange sixth-order basis function
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
            end
 
        case 'h'  % Hermite type FEM element
            if order<3 || order>6
                warning(['The order for Pb & Tb must be an interger more than 2 less than 7 ,',...
                    'already reset to 3']);
                order = 3;
            end
            switch order
                case 3  % Hermite cubic element: a1, da1, a2, da2
                    switch femtype.type
                        case 'pgm' % local basis index for 1st element: 1,2,N+2,N+3
                            
                        case 'ipdg' % local basis index for 1st element: 1,2,3,4
                            
                    end
                case 4  % Hermite fourth-order element: a1, da1, a2, a3, da3
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 5  % Hermite fifth-order element: a1, da1, d^2a1, a2, da2, d^2a2
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 6  % Hermite sixth-order element: a1, da1, d^2a1, a2, a3, da3, d^2a3
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
            end

    end

elseif elemtype == 4  % quadrilateral mesh

    if contains(['l','h','cr','a','b','hct','rt'],type) == 0
        warning(['The type for Pb & Tb must be either l, h, cr, a, b, hct, or rt ',...
            'already reset to l']);
        type = 'l';
    end
    switch type
        case 'l'  % Lagrange type FEM element
            if floor(order)~=order || order<1 || order>7
                warning(['The order for Pb & Tb must be a positive integer less than 8,',...
                    'already reset to 1']);
                order = 1;
            end
            switch order
                case 1 % Lagrange linear basis function
                    switch femtype.type
                        case 'pgm'
                            PT.Pb = mesh.P;
                            PT.Tb = mesh.T;
                        case 'ipdg'
                            
                    end 
                case 2 % Lagrange quadratic basis function
                    switch femtype.type
                        case 'pgm'
                            Pb = mesh.P; Tb = [mesh.T ; mesh.T]; 
                            Tb(4,:) = Tb(2,:); Tb(6,:) = Tb(3,:); 
                            midnode14 = 1/2*(Pb(:,Tb(1,:))+Pb(:,Tb(4,:)));
                            midnode16 = 1/2*(Pb(:,Tb(1,:))+Pb(:,Tb(6,:)));
                            midnode46 = 1/2*(Pb(:,Tb(4,:))+Pb(:,Tb(6,:)));
                            % add midnode between node 1 and node 4 as 2
                            totalnode = size(Pb,2);
                            tmp1 = midnode14(1,:)+1i*midnode14(2,:);
                            [tmp2,~,index] = unique(tmp1,'stable');
                            Pb = [Pb,[real(tmp2);imag(tmp2)]];
                            Tb(2,:) = totalnode + index;
                            % add midnode between node 1 and node 6 as 3
                            totalnode = size(Pb,2);
                            tmp1 = midnode16(1,:)+1i*midnode16(2,:);
                            [tmp2,~,index] = unique(tmp1,'stable');
                            Pb = [Pb,[real(tmp2);imag(tmp2)]];
                            Tb(3,:) = totalnode + index;
                            % add midnode between node 4 and node 6 as 5
                            totalnode = size(Pb,2);
                            tmp1 = midnode46(1,:)+1i*midnode46(2,:);
                            [tmp2,~,index] = unique(tmp1,'stable');
                            Pb = [Pb,[real(tmp2);imag(tmp2)]];
                            Tb(5,:) = totalnode + index;
                            PT.Pb = Pb;
                            PT.Tb = Tb;
                        case 'ipdg'
                            
                    end
                case 3 % Lagrange cubic basis function
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 4 % Lagrange fourth-order basis function
                    switch femtype.type
                        case 'pgm'
                            
                            
                        case 'ipdg'
                            
                    end
                case 5 % Lagrange fifth-order basis function
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 6 % Lagrange sixth-order basis function
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
            end
 
        case 'h'  % Hermite type FEM element
            if order<3 || order>6
                warning(['The order for Pb & Tb must be an interger more than 2 less than 7 ,',...
                    'already reset to 3']);
                order = 3;
            end
            switch order
                case 3  % Hermite cubic element: a1, da1, a2, da2
                    switch femtype.type
                        case 'pgm' % local basis index for 1st element: 1,2,N+2,N+3
                            
                        case 'ipdg' % local basis index for 1st element: 1,2,3,4
                            
                    end
                case 4  % Hermite fourth-order element: a1, da1, a2, a3, da3
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 5  % Hermite fifth-order element: a1, da1, d^2a1, a2, da2, d^2a2
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
                case 6  % Hermite sixth-order element: a1, da1, d^2a1, a2, a3, da3, d^2a3
                    switch femtype.type
                        case 'pgm'
                            
                        case 'ipdg'
                            
                    end
            end

    end

end


end

