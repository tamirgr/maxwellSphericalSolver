sphr =SphereGeometry;        
% background properties
        sphr.ep = 1; %sphr.mu = 1.0;

        % field properties
        sphr.k = 1; %sphr.beta = 0.0;

        % sphrinder properties
        sphr.a = 0.1;
        %sphr.mui = 1;
        %sphr.mu = 2;
        
        % sphrinder coordinates
        indx1 = 10;
        nmodes = 10;
        epiNL1 = zeros(indx1,nmodes);
        epiNL2 = zeros(indx1,nmodes);
         for k = 1:10
             sphr.orders = k-1;
             epiNL1(k,:) = disprootsepi2(sphr, nmodes);
             epiNL2(k,:) = disprootsepi3(sphr, nmodes);
         end
         sphr.orders = 1;
         max(sphrdispTE(sphr,epiNL1(sphr.orders+1,:)))
         max(sphrdispTM(sphr,epiNL2(sphr.orders+1,1:5)))
        % -abs(epiNL2(2,1))^2
         