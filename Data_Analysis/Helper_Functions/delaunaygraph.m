function adjmat = delaunaygraph(xypoints,cutoff,figflag)

% xypoints = coords;
% cutoff = 50;

tri = delaunay(xypoints(:,1),xypoints(:,2));

trilist = tri(:);
trinext = reshape(tri(:,[2 3 1]),size(trilist));
tridist = sqrt(sum((xypoints(trilist,:) - xypoints(trinext,:)).^2,2));

edgeprune = tridist<cutoff;

g = digraph(trilist(edgeprune),trinext(edgeprune),tridist(edgeprune));
A = adjacency(g);%,'weighted');

A = (A + A')/2;




adjmat = logical(A);

%comment out the last line if you don't want ot display the graph
if figflag == 1
    figure
    H = plot(graph(adjmat)); H.XData = xypoints(1:size(adjmat,1),1);H.YData = xypoints(1:size(adjmat,1),2);
end