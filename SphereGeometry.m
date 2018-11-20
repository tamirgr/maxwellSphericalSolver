classdef SphereGeometry
  properties
    % simulation properties
    orders@double;
    ordersN@double;
    
    % background properties
    ep@double scalar complex;

    % field properties
    k@double scalar complex;
    %phi@double scalar complex;
    %psi@double scalar complex;

    % sphere properties
    a@double scalar;
    %epi@double scalar complex;
    
    % Sphere coordinates
    x@double scalar; 
    y@double scalar;
    z@double scalar;

    % sphere coefficients
    AmE; BmE; 
    AmH; BmH; 
  end

  methods
    % index
    function val = n(o)
      val = sqrt(o.ep*o.mu);
    end

    function val = ni(o)
      val = sqrt(o.epi*o.mui);
    end

    % in-plane propagation constant
    function val = kperp(o)
    % val = sqrt([o.k].^2.*[o.ep].*[o.mu] - [o.beta].^2);
      val = sqrtbr([o.k].*sqrt([o.ep].*[o.mu]) - [o.beta], pi) .* sqrtbr([o.k].*sqrt([o.ep].*[o.mu]) + [o.beta], 0);
    end

    function val = kperpi(o)
      val = sqrt([o.k].^2.*[o.epi].*[o.mui] - [o.beta].^2);
    end

    % sphere coordinates polar
    function rc = r(o)
      [~, rc] = cart2pol([o.x], [o.y]);
    end

    function thc = th(o)
      [thc, ~] = cart2pol([o.x], [o.y]);
    end

    % initialize empty coefficients
    function o = initcoeff(o)
      if isempty(o.AmE), o.AmE = zeros(size(o.orders)); end;
      if isempty(o.BmE), o.BmE = zeros(size(o.orders)); end;
      if isempty(o.AmH), o.AmH = zeros(size(o.orders)); end;
      if isempty(o.BmH), o.BmH = zeros(size(o.orders)); end;
    end
  end
end