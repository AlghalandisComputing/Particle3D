% particles3d development bench
% (c) 2020 Alghalandis Computing - https://alghalandis.net
% author: Dr. Younes Fadakar Alghalandis
% last update: 2020-11-27 by C. Masciopinto (CNR)
% to include real case study as "d" file
%

adfnepath = 'ADFNE15';                            								% update this to adfne1.5 path
cd(fileparts(mfilename('fullpath')));                                           % sets the current path
addpath(adfnepath);                                                             % updates Matlab search path
Globals;                                                                        % globals stuff

% ============================================================== reproducibility
rng(123);                                                                       % for reproducibility purpose

% ====================================================================== options
study = 'd';                                                                    % study: r:regular, d:demo
regnn = [21,12];                                                                  % regular case size

% ================================================================== study cases
dem = false;
switch study
    case {'demo','d'}                                                           % demo setup
%%        nodes = [                                                               % nodes (n,3), 3D
%%            0.5	0	10
%%            1	1.5	7
%%            1.5	4	4
%%            2	3	2
%%            8	4	1
%%            7	9.5	0.5
%%            5	9	3
%%            4	7.5	6
%%            3	7	9
%%            5	6	2.5
%%            2.5	6.5	3.5
%%            5	12	2.6
%%            3	12	5
%%            0.5	12	11];
%%        edges = [                                                               % edges (indices to nodes)
%%            1	2
%%            2	3
%%            3	4
%%            4	5
%%            5	6
%%            6	7
%%            7	8
%%            8	9
%%            2	9
%%            9	14
%%            3	11
%%            11	8
%%            8	13
%%            4	10
%%            10	11
%%            10	7
%%            7	12
%%            5	10
%%            6	10];
%%        nodes = Scale(nodes);                                                   % fits to [0,1]
%%
set1 = Field(DFN('dim',3,'n',700,'dir',25,'ddir',-1e9,'minl',0.05,...
'mu',0.1,'maxl',0.5,'bbx',[0,0,0,1,1,1],'dip',45,'ddip',-1e7),'Poly');
set2 = Field(DFN('dim',3,'n',1300,'dir',-15,'ddir',-1e9,'minl',0.05,...
'mu',0.1,'maxl',0.8,'bbx',[0,0,0,1,1,1],'dip',10,'ddip',-1e7,...
'shape','q','q',4),'Poly');
set3 = Field(DFN('dim',3,'n',200,'dir',30,'ddir',-45,'minl',0.05,...
'mu',0.1,'maxl',0.8,'bbx',[0,0,0,1,1,1],'dip',90,'ddip',-1e7,...
'shape','q','q',4),'Poly');
%fnm = Field(DFN('dim',3,'n',200),'Poly');
fnm = [set1;set2;set3];
save('fnm3aab','fnm');
be = {Poly.Left;Poly.Right};
dfn = Pipe(be,fnm);
%%    [eds,edges] = P3FindEdge(dfn,paths);
        dem = true;
    case {'regular','r'}                                                        % regular grid setup
        [nodes,edges] = Grid(regnn(1),regnn(2));
%%end
id = 1:length(nodes);                                                           % id for nodes
of = 1e-6;
be = {Poly.Front-[0,of,0];Poly.Back+[0,of,0]};                                  % boundary entities
cts = Center(be);                                                               % center of boundaries
% =============================================================== pipe structure
nodes = [cts;nodes];                                                            % adds centers to nodes
[~,idx] = sortrows(nodes);                                                      % finds index of updated nodes
ben = [-1,0];                                                                   % boundary labels for nodes
id2 = [ben,id];                                                                 % all node ids
lbl = id2(idx);                                                                 % finds original ids
%if dem                                                                          % demo case
%%edges = [-1,1; 0,12; 0,13; 0,14; edges]+2;                                  % connects center of inlet/outlet
%pn = size(edges,1);                                                         % number of pipes
%typ(1) = 1;                                                                 % 1 pipe, inlet
%typ(2:4) = 2;                                                               % 2..4 pipes, outlet
%else                                                                            % regular case
    edges = [-1,1; 0,length(id); edges]+2;                                      % connects center of inlet/outlet
    pn = size(edges,1);                                                         % number of pipes
    typ = zeros(pn,1);                                                          % pipe types
    typ(1) = 1;                                                                 % 1 pipe, inlet
    typ(2) = 2;                                                                 % 2 pipe, outlet
%end
dfn.pip.Pipe = [nodes(edges(:,1),:),nodes(edges(:,2),:)];                       % building pipe structure...
dfn.pip.n = pn;                                                                 % number of pipes
dfn.pip.Type = typ;                                                             % types
dfn.pip.Label = zeros(pn,1)+1;                                                  % labels
dfn.pip.Depth = ones(pn,1);                                                     % depth for pipes
dfn.pip.Percolate = 1;                                                          % percolation state
dfn.pip.BE = be;                                                                % boundary elements
dfn.dfn.Center = [0,0,0];                                                       % dummy data
end
% ===================================================================== solution
dfn = Solve(Graph(Backbone(dfn)),[0.14,0]);                                        % pressure inlet=1, outlet=0
%Draw(dfn,'what','seniqf'); P3Axes(true); view(45,55); return
   if dem
nodes = Stack({dfn.grh.Node.Location});
edges = Stack({dfn.grh.Edge.Nodes});
ben = [-1,0];
id = 1:length(nodes);
[~,idx] = sortrows(nodes);                                                      % finds index of updated nodes
id2 = [ben,id];                                                                 % all node ids
lbl = id2(idx);
   end
% ====================================================================== options
sni = 12;                                                                        % starting node(s)
npc = 125;                                                                        % number of particles
nts = 80;                                                                       % number of time steps
mxt = 0;                                                                        % <0:min, 0:max, 0<:given
vis = 'a';                                                                      % visual: s:still, a:animate
anmTail = true;                                                                 % to draw tails of particles
anmTitle = true;

% ====================================================================== process
ids = 1:dfn.grh.Size(1);                                                        % counter for nodes
cc = 0.5;
for sn = sni                                                                    % starting node id
    cc = 1.3*cc;
    idx = repmat(ids(lbl == sn),1,npc);                                         % repeats the node id
    
    % ================================================================= pathways
    [eds,paths] = P3FindPath(dfn,ben,lbl,idx);                                    % finds all pathways
    if isempty(eds)                                                             % skips if no pathways
        disp('[!] skipped.');                                                   % ...e.g., at outlet nodes
        continue;
    end
    [ptm,tmc,nod,edg] = P3PathTime(dfn,eds,paths,1.0);                           % computes all arrival times 1.0 is scale factor
    
    % ================================================================= tracking
    if mxt > 0
        tim = linspace(0,mxt,nts);                                              % given time span
    else
        if mxt == 0
            tim = linspace(0,max(ptm),nts);                                     % max time span
        else
            tim = linspace(0,min(ptm),nts);                                     % min time span
        end
    end
    [trk,tci] = P3Track(paths,nod,tmc,tim,lbl);                                   % tracks particles
    
    % =========================================================== visualizations
    b = ([dfn.grh.Edge.Type] == 0);                                             % only inner edges
    P3Plot(nod,edg(b,:),1,1,[0.5,0.5,0.5]);                                     % draws edges and nodes
    b = (lbl > 0);                                                              % only inner nodes
    P3Text(nod(b,:)+[0,0,0.03],lbl(b),4);                                       % draws node labels
    P3Axes(true);
    if ~dem; view(90,90); end
    ntk = length(trk);                                                          % number of tracks
    of = linspace(0,1,npc);                                                     % offset for drawing particles
    of = (of-mean(of))*0.015;
    col = lines(ntk);                                                           % colormap for tracks
    switch vis
        case {'still','s'}                                                      % still image
            for i = 1:ntk
                ci = col(i,:)*cc;
                ci(ci > 1) = 1;
                P3Plot(trk{i}+[0,0,of(i)],[],40,10,ci,0,1);
            end
            if anmSave
                rgb = getframe;
                imwrite(rgb.cdata,'tracking.png');                              % exports as image file, overwrite
            end
        case {'animation','a'}                                                  % animation
            h = [];
            for i = 1:length(idx)                                               % all particles at start
                h = [h,scatter3(0,0,0,150,'filled','MarkerFaceColor',...
                    col(i,:),'MarkerEdgeColor','k')];
            end
            for j = 1:nts                                                       % for all time steps
                for i = 1:ntk                                                   % for all tracks
                    pt = trk{i};                                                % tracks' particle location
                    if j > size(pt,1); continue; end                            % if track already completed
                    if anmTail                                                  % tp draw tailing particles
                        scatter3(pt(j,1),pt(j,2),pt(j,3)+of(i),40,...
                            'filled','MarkerFaceColor',col(i,:));
                    end
                    set(h(i),'XData',pt(j,1),'YData',pt(j,2),'ZData',pt(j,3));  % update particles' location
                end
                if anmTitle; title(['Time: ',num2str(tim(j),'%0.2f')]); end     % title
                drawnow;
            end
    end
end
