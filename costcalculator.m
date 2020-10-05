function [ cost ] = costcalculator (u,n,x,idx)
no_samples = length(x);
dim = size(x,1);

% create voronoi regions
for cin = 1:n
    dum = find(idx==cin);
    voronoiSizes(cin) = length(dum);
    voronois(cin) = mat2cell(x(dum),dim,voronoiSizes(cin));
end

costtemp = 0;
for cin = 1:n
    costi = sum(sqrt((repmat(u(cin), 1, voronoiSizes(cin)) - cell2mat(voronois(cin))).^2))
    costtemp = costi + costtemp;
end
cost = costtemp/no_samples;
end
