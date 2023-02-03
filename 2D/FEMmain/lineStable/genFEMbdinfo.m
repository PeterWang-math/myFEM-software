function trial_test = genFEMbdinfo(bdinfo,type,order,PT,elemtype,femtype)
% bdinfo: boundary information for mesh data
% type: type of FEM basis, the following type is considered:
%       'l' for Lagrange basis funciton,
%       'h' for Hermite basis function,
%       'cr' for Crouzeix-Raviart element,
%       'a' for Argyris element,
%       'b' for Bell element,
%       'hct' for Hsich-Clough-Tocher element,
%       'rt' for Raviart-Thomas element
% order: a positive integer representing the order of basis function, 
%        less than 8 
% PT: PT.Pb & PT.Tb are FEM information matrix
% elemtype: 3: triangular mesh; 4: quadrilateral mesh
% femtype: 'pgm': Petrov-Galerkin method; 
%          'ipdg': internal penalty discontinous Galerkin method;

if elemtype == 3  % triangular mesh
    switch femtype.type
        case 'pgm' % Petrov-Galerkin method
            switch type
                case 'l'
                    switch order
                        case 1  % Lagrange linear element
                            trial_test.edges = bdinfo.edges;
                            trial_test.nodes = bdinfo.nodes;
                        case 2  % Lagrange quadratic element
                            % find global number of midnode
                            newnode = zeros(2,size(bdinfo.edges,2));
                            newnode(1,:) = bdinfo.edges(1,:);
                            midnode = ( PT.Pb(:,bdinfo.edges(3,:)) + PT.Pb(:,bdinfo.edges(4,:)) )/2;
                            [~,index] = ismember(midnode(1,:)+1i*midnode(2,:),PT.Pb(1,:)+1i*PT.Pb(2,:));
                            newnode(2,:) = index;
                            trial_test.nodes = [bdinfo.nodes, newnode];
                            % get new edge
                            newedge = [bdinfo.edges, bdinfo.edges];
                            newedge(3,:) = [bdinfo.edges(3,:), index];
                            newedge(4,:) = [index, bdinfo.edges(4,:)];
                            trial_test.edges = newedge;
                    end

                case 'h'
                    switch order
                        case 3

                        case 4

                    end
            end

        case 'ipdg' % internal penalty discontinous Galerkin method
            switch type
                case 'l'
                    switch order
                        case 1

                        case 2

                    end

                case 'h'
                    switch order
                        case 3

                        case 4

                    end
            end
    end
elseif elemtype == 4  % quadrilateral mesh
    switch femtype.type
        case 'pgm' % Petrov-Galerkin method
            switch type
                case 'l'
                    switch order
                        case 1  % Lagrange linear element
                            trial_test.edges = bdinfo.edges;
                            trial_test.nodes = bdinfo.nodes;
                        case 2  

                    end

                case 'h'
                    switch order
                        case 3

                        case 4

                    end
            end

        case 'ipdg' % internal penalty discontinous Galerkin method
            switch type
                case 'l'
                    switch order
                        case 1

                        case 2

                    end

                case 'h'
                    switch order
                        case 3

                        case 4

                    end
            end
    end
end


end

