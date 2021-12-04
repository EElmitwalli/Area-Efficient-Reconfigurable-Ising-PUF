function [APUFr] = APUFResponseTable(APUFw,maxconn)


%Generate binary truthtable 
nRows = 2^maxconn;
N=maxconn;
L = 2^N;
TeS = zeros(L,N);
for i=1:N
   temp = [zeros(L/2^i,1); ones(L/2^i,1)];
   TeS(:,i) = repmat(temp,2^(i-1),1);
end

%Transform to feature vector
APhi = Transform(TeS, nRows, maxconn);

%Calculate PUF response
temp=APhi*APUFw';
APUFr=temp >= 0;
