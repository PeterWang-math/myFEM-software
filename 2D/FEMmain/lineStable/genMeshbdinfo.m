function newbdinfo = genMeshbdinfo(varargin)
% generate boundary edges and nodes information for input mesh data
% the input elem & bdFlag is from iFEM software package
% newbdinfo = genMeshbdinfo(elem,bdFlag)
% newbdinfo = genMeshbdinfo(elem,bdFlag,bdinfo)

if nargin == 2
    elem = varargin{1};
    bdFlag = varargin{2};
elseif nargin == 3 && isfield(varargin{3},'function')  % store bdinfo.function
    elem = varargin{1};
    bdFlag = varargin{2};
    bdinfo = varargin{3};
    newbdinfo.function = bdinfo.function;
end

switch size(elem,2)
    case 3  % triangular mesh
        % generate boundary edges' information
        allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        bdEdge = allEdge((bdFlag(:) ~= 0),:); % the index of all boundary edges
        Dirichlet = allEdge((bdFlag(:) == 1),:); % Dirichlet b.c.
        dlen = size(Dirichlet,1);
        Neumann = allEdge((bdFlag(:) == 2),:); % Neumann b.c.
        nlen = size(Neumann,1);
        Robin = allEdge((bdFlag(:) == 3),:); % Robin b.c.
        rlen = size(Robin,1);
        edge2bd = [Dirichlet,ones(dlen,1); Neumann,2*ones(nlen,1); Robin,3*ones(rlen,1)];
        edgelen = size(bdEdge,1);
        elemlen = size(elem,1);
        edge2elem = [elem(:,[2,3]),(1:elemlen)'; ...
            elem(:,[3,1]),(1:elemlen)'; ...
            elem(:,[1,2]),(1:elemlen)'];
        [lia, locb] = ismember(allEdge,edge2bd(:,1:2),'row'); locb(locb==0)=[];
        edgeFlag = edge2bd(:,3);

        bdinfo.edges = zeros(4,edgelen);
        bdinfo.edges(1,:) = edgeFlag(locb);
        bdinfo.edges(2,:) = edge2elem(lia==1,3);
        bdinfo.edges(3:4,:) = bdEdge';

        % generate boundary nodes' information
        [bdNode,~,~] = findboundary(elem);
        nodelen = length(bdNode);
        dnode = unique(Dirichlet); dlen = size(dnode,1);
        nnode = unique(Neumann); nlen = size(nnode,1);
        rnode = unique(Robin); rlen = size(rnode,1);
        node2bd = [dnode, ones(dlen,1); nnode, 2*ones(nlen,1); rnode, 3*ones(rlen,1)];
        [~, locb] = ismember(bdNode,node2bd(:,1),'row');
        nodeFlag = node2bd(:,2);

        bdinfo.nodes = zeros(2,nodelen);
        bdinfo.nodes(1,:) = nodeFlag(locb)';
        bdinfo.nodes(2,:) = bdNode';

    case 4  % quadrilateral mesh
        % for quadrilateral element, the index number must be counterclockwise
        %t = elem(:,2); elem(:,2) = elem(:,4); elem(:,4) = t;

        % generate boundary edges' information
        allEdge = [elem(:,[1,2]); elem(:,[2,3]); elem(:,[3,4]);elem(:,[4,1])];
        bdEdge = allEdge((bdFlag(:) ~= 0),:); % the index of all boundary edges
        Dirichlet = allEdge((bdFlag(:) == 1),:); % Dirichlet b.c.
        dlen = size(Dirichlet,1);
        Neumann = allEdge((bdFlag(:) == 2),:); % Neumann b.c.
        nlen = size(Neumann,1);
        Robin = allEdge((bdFlag(:) == 3),:); % Robin b.c.
        rlen = size(Robin,1);
        edge2bd = [Dirichlet,ones(dlen,1); Neumann,2*ones(nlen,1); Robin,3*ones(rlen,1)];
        edgelen = size(bdEdge,1);
        elemlen = size(elem,1);
        edge2elem = [elem(:,[1,2]),(1:elemlen)'; ...
            elem(:,[2,3]),(1:elemlen)'; ...
            elem(:,[3,4]),(1:elemlen)'; ...
            elem(:,[4,1]),(1:elemlen)'];
        [lia, locb] = ismember(allEdge,edge2bd(:,1:2),'row'); locb(locb==0)=[];
        edgeFlag = edge2bd(:,3);

        bdinfo.edges = zeros(4,edgelen);
        bdinfo.edges(1,:) = edgeFlag(locb);
        bdinfo.edges(2,:) = edge2elem(lia==1,3);
        bdinfo.edges(3:4,:) = bdEdge';

        % generate boundary nodes' information
        [bdNode,~,~] = findquadboundary(elem);
        nodelen = length(bdNode);
        dnode = unique(Dirichlet); dlen = size(dnode,1);
        nnode = unique(Neumann); nlen = size(nnode,1);
        rnode = unique(Robin); rlen = size(rnode,1);
        node2bd = [dnode, ones(dlen,1); nnode, 2*ones(nlen,1); rnode, 3*ones(rlen,1)];
        [~, locb] = ismember(bdNode,node2bd(:,1),'row');
        nodeFlag = node2bd(:,2);

        bdinfo.nodes = zeros(2,nodelen);
        bdinfo.nodes(1,:) = nodeFlag(locb)';
        bdinfo.nodes(2,:) = bdNode';

end

newbdinfo.edges = bdinfo.edges;
newbdinfo.nodes = bdinfo.nodes;

end

