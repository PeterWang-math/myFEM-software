function trial_test = genTribasis2(type,order)
% generate trial_test basis function space infomation for triangular element
% trial_test space is considered seperately for PGM
%
% The following input infomation is needed:
% type: type of FEM basis, the following type is considered:
%       'l' for Lagrange basis funciton,
%       'h' for Hermite basis function,
%       'cr' for Crouzeix-Raviart element,
%       'a' for Argyris element,
%       'b' for Bell element,
%       'hct' for Hsich-Clough-Tocher element,
%       'rt' for Raviart-Thomas element
% order: a positive integer representing the order of basis function, less than 8
%
% The output information includes:
% trial_test.reffun: Nlb x (num_of_derivatives+1) x (num_of_derivatives+1) 
%                     cell consisting of all the local basis functions on 
%                     reference element as well as its derivatives, the 
%                     {k,s+1,s+t+1} stores d^{s+t}phi_k/dx^{s}dy^{t}
% e.g.: s+t=3, yyy,xyy,xxy,xxx
% trial_test.Nb/Nlb: FEM nodes information

% linear nodal basis
%l1 = @(x,y) 1-x-y;
%l2 = @(x,y) x;
%l3 = @(x,y) y;

trial_test.type = type;
trial_test.order = order;

%% generate trial_test basis function
% the basis function is defined on the reference triangle element
if contains(['l','h','cr','a','b','hct','rt'],trial_test.type) == 0
    warning(['The trial_test.type must be either l, h, cr, a, b, hct, or rt ',...
        'already reset to l']);
    trial_test.type = 'l';
end
%% vectorization programming: every constant C in fclass is modified as C+x-x
switch trial_test.type
    case 'l'  % Lagrange type FEM element
        ord = trial_test.order;
        if floor(ord)~=ord || ord<1 || ord>7
            warning(['The trial_test.order must be a positive integer less than 8,',...
                'already reset to 1']);
            trial_test.order = 1;
        end
        switch trial_test.order
            case 1 % Lagrange linear basis function
                fclass = cell(3,8,8);
                %  fclass{1,1,1} = @(x,y) l1(x,y);
                %  fclass{2,1,1} = @(x,y) l2(x,y);
                %  fclass{3,1,1} = @(x,y) l3(x,y);
                fclass{1,1,1} = @(x,y) 1-x-y;
                fclass{2,1,1} = @(x,y) x;
                fclass{3,1,1} = @(x,y) y;
                fclass{1,1,2} = @(x,y) -1+x-x;
                fclass{2,1,2} = @(x,y) 0+x-x;
                fclass{3,1,2} = @(x,y) 1+x-x;
                fclass{1,2,2} = @(x,y) -1+x-x;
                fclass{2,2,2} = @(x,y) 1+x-x;
                fclass{3,2,2} = @(x,y) 0+x-x;
                for i=3:8
                    for j=1:3
                        for k=1:8
                            if k<=i
                                fclass{j,k,i} = @(x,y) 0+x-x;
                            end
                        end
                    end
                end
                trial_test.reffun = fclass;
            case 2 % Lagrange quadratic basis function
                fclass = cell(6,8,8);
                fclass{1,1,1} = @(x,y) (-1+x+y).*(-1+2*x+2*y);
                fclass{2,1,1} = @(x,y) -4*x.*(-1+x+y);
                fclass{3,1,1} = @(x,y) -4*y.*(-1+x+y);
                fclass{4,1,1} = @(x,y) x.*(-1+2*x);
                fclass{5,1,1} = @(x,y) 4*x.*y;
                fclass{6,1,1} = @(x,y) y.*(-1+2*y);
                % uy
                fclass{1,1,2} = @(x,y) -3+4*x+4*y;
                fclass{2,1,2} = @(x,y) -4*x;
                fclass{3,1,2} = @(x,y) -4.*(-1+x+2*y);
                fclass{4,1,2} = @(x,y) 0+x-x;
                fclass{5,1,2} = @(x,y) 4*x;
                fclass{6,1,2} = @(x,y) -1+4*y;
                % ux
                fclass{1,2,2} = @(x,y) -3+4*x+4*y;
                fclass{2,2,2} = @(x,y) -4*(-1+2*x+y);
                fclass{3,2,2} = @(x,y) -4*y;
                fclass{4,2,2} = @(x,y) -1+4*x;
                fclass{5,2,2} = @(x,y) 4*y;
                fclass{6,2,2} = @(x,y) 0+x-x;
                % uyy
                fclass{1,1,3} = @(x,y) 4+x-x;
                fclass{2,1,3} = @(x,y) 0+x-x;
                fclass{3,1,3} = @(x,y) -8+x-x;
                fclass{4,1,3} = @(x,y) 0+x-x;
                fclass{5,1,3} = @(x,y) 0+x-x;
                fclass{6,1,3} = @(x,y) 4+x-x;
                % uxy
                fclass{1,2,3} = @(x,y) 4+x-x;
                fclass{2,2,3} = @(x,y) -4+x-x;
                fclass{3,2,3} = @(x,y) -4+x-x;
                fclass{4,2,3} = @(x,y) 0+x-x;
                fclass{5,2,3} = @(x,y) 4+x-x;
                fclass{6,2,3} = @(x,y) 0+x-x;
                % uxx
                fclass{1,3,3} = @(x,y) 4+x-x;
                fclass{2,3,3} = @(x,y) -8+x-x;
                fclass{3,3,3} = @(x,y) 0+x-x;
                fclass{4,3,3} = @(x,y) 4+x-x;
                fclass{5,3,3} = @(x,y) 0+x-x;
                fclass{6,3,3} = @(x,y) 0+x-x;
                for i=4:8
                    for j=1:6
                        for k=1:8
                            if k<=i
                                fclass{j,k,i} = @(x,y) 0+x-x;
                            end
                        end
                    end
                end
                trial_test.reffun = fclass;
            case 3 % Lagrange cubic basis function
                fclass = cell(7,4);

                trial_test.reffun = fclass;
            case 4 % Lagrange fourth-order basis function
                fclass = cell(7,5);

                trial_test.reffun = fclass;
            case 5 % Lagrange fifth-order basis function
                fclass = cell(7,6);

                trial_test.reffun = fclass;
            case 6 % Lagrange sixth-order basis function
                fclass = cell(7,7);

                trial_test.reffun = fclass;
        end


    case 'h'   % Hermite type FEM element
        ord = trial_test.order;
        if ord<3 || ord>6
            warning(['The trial_test.order must be an interger more than 2 less than 7 ,',...
                'already reset to 3']);
            trial_test.order = 3;
        end
        switch trial_test.order
            case 3  % Hermite cubic element: a1, da1, a2, da2
                fclass = cell(7,4);

                trial_test.reffun = fclass;
            case 4  % Hermite fourth-order element: a1, da1, a2, a3, da3
                fclass = cell(7,5);

                trial_test.reffun = fclass;
            case 5  % Hermite fifth-order element: a1, da1, d^2a1, a2, da2, d^2a2
                fclass = cell(7,6);

                trial_test.reffun = fclass;
            case 6  % Hermite sixth-order element: a1, da1, d^2a1, a2, a3, da3, d^2a3
                fclass = cell(7,7);

                trial_test.reffun = fclass;
        end

end


end
