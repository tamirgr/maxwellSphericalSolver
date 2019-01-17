function [valx, valy, valz] = greentensorint(rs, ths, zs, xd, yd, zd, k, inx, iny, inz)

% integral of volume source with Helmholtz Green's tensor kernel
% cylindrical domain: midpoint rule
% non-coincident source and detector points
% source: cylindrical, detector: cartesian

% collect inputs
in{1} = inx; in{2} = iny; in{3} = inz;

% step sizes
dr = (rs(end)-rs(1))./(length(rs)-1);
dth = 2.*pi./length(ths);
dz = (zs(end)-zs(1))./(length(zs)-1);
dV = dr.*dth.*dz;

% construct arrays of quadrature (source) points
[r, th, z] = ndgrid(rs, ths, zs);

% to obtain cartesian distance
[x, y] = pol2cart(th, r);

% initialize results
[val{1:3}] = deal(complex(zeros(size(xd))));

% iterate over detector points
for it = 1:numel(xd)
  % components of distance separation
  sep{1} = x-xd(it);
  sep{2} = y-yd(it);
  sep{3} = z-zd(it);
  
  % distance separation
  dist = sqrt(sep{1}.^2 + sep{2}.^2 + sep{3}.^2);
  
  % separation unit vectors
  uvec{1} = sep{1}./dist;
  uvec{2} = sep{2}./dist;
  uvec{3} = sep{3}./dist;
  
  % Green's scalar w/ Jacobian
  gscalar = r .* exp(1i.*k.*dist)./dist;
  
  % Green's tensor polynomial terms
  inv = 1./k./dist;
  invs = inv.^2;
    
  % midpoint rule integration
  % iterate over 9 tensorial components
  for s = 1:3
    for d = 1:3
      val{d}(it) = val{d}(it) + sum(sum(sum((3.*invs-3i.*inv-1).*uvec{s}.*uvec{d}.*gscalar.*in{s}))).*dV;
%     val{d}(it) = val{d}(it) + sum(sum(sum(-uvec{s}.*uvec{d}.*gscalar.*in{s}))).*dV;
      if s == d
%       val{d}(it) = val{d}(it) + sum(sum(sum(gscalar.*in{s}))).*dV;
        val{d}(it) = val{d}(it) + sum(sum(sum((1+1i.*inv-invs).*gscalar.*in{s}))).*dV;
      end
    end
  end
end

% distribute outputs
valx = val{1}; valy = val{2}; valz = val{3};
