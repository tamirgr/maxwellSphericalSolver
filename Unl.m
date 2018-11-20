function [res] = Unl(sphr,r)
    ext = r > sphr.a;
    res = ((sphr.ep - sphr.epi)./shpr.ep).*ext;
end