%% Vessel set-up
numNodes = size(coord,2); % Total number of nodes in the model

numElem = size(connect,2); % Total number of elements in the model

numDofPerNode = size(coord,1); % Number of velocity components

numNodesPerElem = size(connect,1); % Number of nodes per element

numDofPerElem = numNodesPerElem * numDofPerNode; 
                                 % Number of DOF associated with an element
ModelSize = (numNodes*numDofPerNode + 3*(N+1)^2);

%% Set the nodal boundary conditions
% The nodal values of the zero-velocity field on the sidewall surface are
% assigned with the Dirichlet nodes, which are the nodes on the sidewall
% including the edge nodes.
% Similarly, the nodal values of the inlet pressure and zero-pressure
% fields are assigned with the Neumann nodes, which are the nodes on the
% inlet and the outlet surfaces excluding the edge nodes.
% Therefore, the unknown fields consist of the nodal values of the inlet
% velocity, the outlet velocity and the wall traction fields.

% bcs == 1; Dirichlet B.C; velocity is known
% bcs == 0; Neumann B.C; traction is known
bcs = zeros(size(coord));
bcs(:) = 1;
bcs(:,setdiff(inletnode,wallnode)) = 0;
bcs(:,setdiff(outletnode,wallnode)) = 0;

DirichletNode = wallnode;
NeumannNode = union(setdiff(inletnode,wallnode), ...
                    setdiff(outletnode,wallnode));

DirichletDofs = find(bcs == 1); % Dofs of the wall nodes
NeumannDofs = find(bcs == 0);   % Dofs of the inlet/outlet nodes 
                                % excluding dofs related to the edge nodes.

%% Set element boundary conditions
bcselem = zeros(1, size(connect,2));
bcselem(1,wallelem) = 1; % Dirichlet B.C; velocity is known

DirichletElem = wallelem; % = find(bcselem);
NeumannElem = union(inletelem,outletelem); % = find(~bcselem);

%% Set prescribed (inlet) traction
Telem = zeros(numDofPerElem, numElem);
%% Prescribe the inlet pressure on the inlet element nodes in x-direction
Telem(1:numDofPerNode:numDofPerElem,inletelem) = InletPressure; 

%% Set prescribed velocity (all zero for now)
unodal = zeros(numDofPerNode, numNodes);

%% Array to store element DOF numbers
elemDofNum = zeros(numDofPerElem, numElem);
for m = 1:numElem
    for n = 1:numNodesPerElem
        elemDofNum((n-1)*numDofPerNode+(1:numDofPerNode),m) = ...
                          (connect(n,m)-1)*numDofPerNode+(1:numDofPerNode);
    end
end