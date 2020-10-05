function [ cost ] = ObjectiveFunction ( u,n,r,x,h,delta,P,b,c )
no_samples = length(x);
dim = size(x,1);
for cin = 1:n
    csDists(:,cin) = sqrt(sum((repmat(u(:,cin), 1, no_samples) - x).^2, 1));
end
% nearest neighbors of each sample
[dist_to_nns, nnss] = min(csDists,[],2);
if n==1
    nns=1;
else
    nns=nnss;
end

% create voronoi regions
for cin = 1:n
    dum = find(nns==cin);
    voronoiSizes(cin) = length(dum);
    voronois(cin) = mat2cell(x(:,dum),dim,voronoiSizes(cin));
end

% Calculate the cost
costtemp = 0;
for cin = 1:n
    PLOS = 1./(1+c*exp(-b*(atan(h./sqrt(sum((repmat(u(:,cin), 1, voronoiSizes(cin)) - cell2mat(voronois(cin))).^2, 1)))-c))) ;
    costLOS = sum(log2(1+(P./(sum((repmat(u(:,cin), 1, voronoiSizes(cin)) - cell2mat(voronois(cin))).^2, 1)+h^2).^(r/2))).*PLOS);
    costNLOS = sum(log2(1+(P*delta./(sum((repmat(u(:,cin), 1, voronoiSizes(cin)) - cell2mat(voronois(cin))).^2, 1)+h^2).^(r/2))).*(1-PLOS));
    costtemp = costLOS + costNLOS + costtemp;
end
cost = -costtemp/no_samples;
end
