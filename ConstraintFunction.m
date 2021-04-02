clc;
clear all;
global E
global FS
global Rho
global YS
global Allowable_displacement
global Coordinates
global ElementCon
global NumberElement
global NumberNode
global Fg
global U_glob

E = 6900;
FS = 1.5;
Rho = 2.77e-5;
YS = 172;
Allowable_displacement = 50.8;
Coordinates = [0 0; 9144 0;18288 0; 18288 9144; 9144 9144;0 9144];
ElementCon = [1 2; 2 3; 3 4; 4 5; 5 6;2 5; 2 6; 1 5; 3 5; 2 4];
NumberNode = size(Coordinates,1);
NumberElement = size(ElementCon,1);
Kg = zeros(2*NumberNode);
A = 3.14*ones(1,NumberElement);
 
ForceBC = input('Enter the total number of nodes on which force acts\n');
    Fg = zeros(2*NumberNode,1);
    for i = 1:ForceBC
        fprintf('Enter the Node number on which force %d acts\n',i);
        NodeFBC = input('');
        Fnode = input('Enter the x, y component of force in pounds on the node in [Fx;Fy] format\n');
        Fg((NodeFBC*2)-1:(NodeFBC*2),1) = Fg((NodeFBC*2)-1:(NodeFBC*2),1) + Fnode;
    end
U_glob = ones(2*NumberNode,1);
DisplacementBC = input('Enter the total number of nodes on which displacement constraints are applied\n');
for i = 1:DisplacementBC
    NodeDBC = input('Enter the Node ID  = ');
    U_glob(2*NodeDBC-1:2*NodeDBC) = 0;
end

b = find(U_glob == 0);
Fg(b,:) = [];

Fitness = @ObjectiveFunction;
Constraint = @ConstraintFunction;
nvars = 10;
LB= 0.1*A;
UB = 2*A;
% opts=gaoptimset('PopulationSize',20,'Generations',200,'StallGenLimit',200,'SelectionFcn', @selectionroulette,'CrossoverFcn',@crossovertwopoint,'Display', 'off','PlotFcns', {@gaplotbestf @gaplotbestindiv});
opts = optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotbestindiv});
disp('Working on optimizing');
[x,fval] = ga(Fitness,nvars,[],[],[],[],LB,UB,Constraint,opts);

function [constraint_ieq,constraint_eq] = ConstraintFunction(A)
    global E
    global FS
    global Rho
    global YS
    global Allowable_displacement
    global Coordinates
    global ElementCon
    global NumberElement
    global NumberNode
    global Fg
    global U_glob
    Kg = zeros(2*NumberNode);
    for i = 1:NumberElement
        n1 = ElementCon(i,1);
        n2 = ElementCon(i,2); 
        x1 = Coordinates(n1,1);
        y1 = Coordinates(n1,2);
        x2 = Coordinates(n2,1);
        y2 = Coordinates(n2,2);

        L = sqrt((x2-x1)^2+(y2-y1)^2);
        l = (x2-x1)/L;
        m = (y2-y1)/L;

        Kl = (E*A(i)/L)*[l^2 l*m -l^2 -l*m;
                        l*m m^2 -l*m -m^2;
                        -l^2 -l*m l^2 l*m;
                        -l*m -m^2 l*m m^2];

        k1 = 2*n1-1; k2 = 2*n1;
        k3 = 2*n2-1; k4 = 2*n2;

        Kg(k1:k2,k1:k2) = Kg(k1:k2,k1:k2) + Kl(1:2,1:2);
        Kg(k1:k2,k3:k4) = Kg(k1:k2,k3:k4) + Kl(1:2,3:4);
        Kg(k3:k4,k1:k2) = Kg(k3:k4,k1:k2) + Kl(3:4,1:2);
        Kg(k3:k4,k3:k4) = Kg(k3:k4,k3:k4) + Kl(3:4,3:4);
    end
    b = find(U_glob == 0);
    Kg(b,:) = [];
    Kg(:,b) = [];
%     Fg(b,:) = [];
    bb = find(U_glob == 1);
    Q = Kg\Fg ;
    Q_disp = zeros(2*NumberNode,1);
    for i=1:length(bb)
      Q_disp(bb(i),1) = Q(i,1);
    end
    for i=1:NumberElement
        L(i)=sqrt((Coordinates(ElementCon(i,2),1)-Coordinates(ElementCon(i,1),1))^2+(Coordinates(ElementCon(i,2),2)-Coordinates(ElementCon(i,1),2))^2);
        l(i)=(Coordinates(ElementCon(i,2),1)-Coordinates(ElementCon(i,1),1))/L(i);
        m(i)=(Coordinates(ElementCon(i,2),2)-Coordinates(ElementCon(i,1),2))/L(i);
    end
    Sigma = zeros(1,NumberElement);
    for k = 1:NumberElement
        Sigma(k) = (E/L(k))*[-1 1]*[l(k), m(k),0, 0;0, 0 , l(k), m(k)] * [Q_disp(((ElementCon(k,1)*2)-1):ElementCon(k,1)*2,1);Q_disp(((ElementCon(k,2)*2)-1):ElementCon(k,2)*2,1)];
    end
    stress = transpose(Sigma);
    max_displacement = max(Q_disp);
    max_stress = max(stress);
    cstress = max_stress - YS;
    cstrain = max_displacement - Allowable_displacement;
    cstress_2 = abs(min(stress)) - YS;
    constraint_ieq = [cstress,cstress_2];
    constraint_eq = [];
end