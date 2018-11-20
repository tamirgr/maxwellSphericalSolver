function [E, H] = fieldpt2(x, y, z, sphr, varargin)

% Obtains cylindrical or Cartesian fields components at specified point

% plot points
%x(x==0 & y==0) = x(x==0 & y==0) + eps; % coordinate singularity
[th, phi, r] = cart2sph(x-sphr.x, y-sphr.y,z-sphr.z);

% test if plot points exterior to spherical
ext = r > sphr.a;

% fill in empty fields
%c = sphr.initcoeff;
%for every value of n,l,m there is a different electrical field
[E, H] = deal(zeros(size(r)));
for n = 1:length(sphr.ordersN)
%    for l = 1:length(sphr.orders);
%   bessj = besselj(sphr.orders(s),sphr.kperp.*r(ext));
%   bessh = besselh1(sphr.orders(s),sphr.kperp.*r(ext));
%   bessjd = besseljd(sphr.orders(s),sphr.kperp.*r(ext));
%   besshd = besselh1d(sphr.orders(s),sphr.kperp.*r(ext));
%   bessji = besselj(sphr.orders(s),sphr.kperpi.*r(~ext));
%   besshi = besselh1(sphr.orders(s),sphr.kperpi.*r(~ext));
%   bessjid = besseljd(sphr.orders(s),sphr.kperpi.*r(~ext));
%   besshid = besselh1d(sphr.orders(s),sphr.kperpi.*r(~ext));
% 
%   expk = exp(i.*sphr.orders(s).*th(ext));
%   expki = exp(i.*sphr.orders(s).*th(~ext));
  
    xlm = Xlm(n,th,phi);
    tmp = shiftdim(repmat(xlm,[ones(1,length(size(xlm))),length(r)]),length(size(xlm)));
    E = E + bsxfun(@times,tmp,transpose(Fnl(n,r,sphr,'e')));
%  Et = Et + 1i./sphr.kperp.^2.*(sphr.beta./r(ext).*i.*sphr.orders(s).*(sphr.AmE(s).*bessj+sphr.BmE(s).*bessh)-sphr.k.*sphr.mu.*sphr.kperp.*(sphr.AmH(s).*bessjd+sphr.BmH(s).*besshd)).*expk;
%  Er = Er + i./sphr.kperp.^2.*(sphr.beta.*sphr.kperp.*(sphr.AmE(s).*bessjd+sphr.BmE(s).*besshd)+sphr.k.*sphr.mu./r(ext).*i.*sphr.orders(s).*(sphr.AmH(s).*bessj+sphr.BmH(s).*bessh)).*expk;

%   Ep(~ext) = Ep(~ext) + (sphr.CmE(s).*bessji+sphr.DmE(s).*bessh  i).*expki;
%   Et(~ext) = Et(~ext) + i./sphr.kperpi.^2.*(sphr.beta./r(~ext).*i.*sphr.orders(s).*(sphr.CmE(s).*bessji+sphr.DmE(s).*besshi) - sphr.k.*sphr.mui.*sphr.kperpi.*(sphr.CmH(s).*bessjid+sphr.DmH(s).*besshid)).*expki;
%   Er(~ext) = Er(~ext) + i./sphr.kperpi.^2.*(sphr.beta.*sphr.kperpi.*(sphr.CmE(s).*bessjid+sphr.DmE(s).*besshid)+sphr.k.*sphr.mui./r(~ext).*i.*sphr.orders(s).*(sphr.CmH(s).*bessji+sphr.DmH(s).*besshi)).*expki;

%    H = H + 1i/k0./bsxfun(@times,(1-Unl(sphr,r)),curlFXlm(Fnl(n,r,sphr,'m'),xlm));
%  Ht = Ht + i./sphr.kperp.^2.*(sphr.beta./r(ext).*i.*sphr.orders(s).*(sphr.AmH(s).*bessj+sphr.BmH(s).*bessh)+sphr.k.*sphr.ep.*sphr.kperp.*(sphr.AmE(s).*bessjd+sphr.BmE(s).*besshd)).*expk;
%  Hr = Hr + i./sphr.kperp.^2.*(sphr.beta.*sphr.kperp.*(sphr.AmH(s).*bessjd+sphr.BmH(s).*besshd)-sphr.k.*sphr.ep./r(ext).*i.*sphr.orders(s).*(sphr.AmE(s).*bessj+sphr.BmE(s).*bessh)).*expk;
%   inside the sphere
%   Hp(~ext) = Hp(~ext) + (sphr.CmH(s).*bessji+sphr.DmH(s).*besshi).*expki;
%   Ht(~ext) = Ht(~ext) + i./sphr.kperpi.^2.*(sphr.beta./r(~ext).*i.*sphr.orders(s).*(sphr.CmH(s).*bessji+sphr.DmH(s).*besshi) + sphr.k.*sphr.epi.*sphr.kperpi.*(sphr.CmE(s).*bessjid+sphr.DmE(s).*besshid)).*expki;
%   Hr(~ext) = Hr(~ext) + i./sphr.kperpi.^2.*(sphr.beta.*sphr.kperpi.*(sphr.CmH(s).*bessjid+sphr.DmH(s).*besshid)-sphr.k.*sphr.epi./r(~ext).*i.*sphr.orders(s).*(sphr.CmE(s).*bessji+sphr.DmE(s).*besshi)).*expki;
    end
end

% Cartesian fields
if nargin > 3 & varargin{1} == 'cart'
  % conversion
  Ex = sin(th).*(cos(phi).*Er - sin(phi).*Et) + cos(phi).*Ep;
  Ey = sin(th).*(sin(phi).*Er + cos(phi).*Et) + sin(phi).*Ep;
  Ez = cos(th).*Er + cos(th).*Et;
  % reassign outputs
  Er = Ex; Et = Ey;
end
