clc,clear
tic

isingSize=4;  %number of APUFs (x) in a side length (i.e. x*x-Ising-PUF)
maxconn=4; %number of connections per node (4-8)

annealing_steps=1;  %annealing steps annealing_steps

%APUF parameters
mu=0.1;
sigma=1;

%number of challenges to generate
CRPCount=1000000;

%reproducable randomness
rseed=12; % RNG seed here
rng(rseed,'combRecursive')

%preallocate arrays
PUFResponse=nan;
challengebitsPhi=zeros(CRPCount,isingSize^2*(maxconn+1));

%generate APUF weights
for i=1:isingSize^2
	APUFw(i,:) = normrnd(mu,sigma,1,maxconn+1);
end

%Calculate APUFs CRPs
APUF=APUFResponseTable(APUFw,maxconn);

%APUF binary conversion variable (saves time to calculate it here)
for i=1:maxconn
    maxconnConVar(i)= 2^(maxconn-i);
end

parfor TrainingSet=1:CRPCount
    %preallocate arrays
    neigh=zeros(1,maxconn);
    conn=zeros(1,maxconn);
    neighbourMatrix=zeros(isingSize^2,maxconn);
    
    isingMachine = randi([0 1], isingSize,isingSize); %random Spin Challenge
    connectionMatrix = randi([0 1], isingSize^2,maxconn); %random Connection Challenge
    
    %Apply annealing
    for steps=1:annealing_steps
        isingState= isingMachine;
        
        %update each bit according to challenge from neighbours
        for m=1:isingSize
            for n=1:isingSize
                
                linIndex=n+isingSize*(m-1);
                
                %get neighbours for each APUF
                neigh(1)=isingMachine(mod(m-2,isingSize)+1,n); %top neighbour
                neigh(2)=isingMachine(m,mod(n,isingSize)+1); %right neighbour
                neigh(3)=isingMachine(mod(m,isingSize)+1,n); %bottom neighbour
                neigh(4)=isingMachine(m,mod(n-2,isingSize)+1); %left neighbour
                
                %get connection for each APUF
                conn=connectionMatrix(linIndex,1:maxconn);
                
                %extra neighbours addition
                if maxconn > 4 %fifth connection
                    neigh(5)=isingMachine(mod(m-2,isingSize)+1,mod(n,isingSize)+1); %top right neighbour
                    
                    if maxconn > 5 %sixth connection
                        neigh(6)=isingMachine(mod(m,isingSize)+1,mod(n,isingSize)+1); %bot right neighbour
                        
                        if maxconn > 6 %seventh connection
                            neigh(7)=isingMachine(mod(m,isingSize)+1,mod(n-2,isingSize)+1); %bot left neighbour
                            
                            if maxconn > 7 %eightth connection
                                neigh(8)=isingMachine(mod(m-2,isingSize)+1,mod(n-2,isingSize)+1); %top left neighbour
                            end
                        end
                    end
                end
                
                %Check if any XOR gates are active
                if sum(conn)~=0
                    neigh=mod(neigh+conn,2);
                end
                
                %Save challenge sequence for each APUF node in first step
                if steps==1
                    neighbourMatrix(linIndex,1:maxconn)=neigh;
                end
                
                %Update next state for this bit
                isingState(m,n)=APUF(sum(neigh.*maxconnConVar)+1,linIndex);
            end
        end
        isingMachine=isingState;
    end
    
    % XOR after # annealing steps to get IsingPUF response
    PUFResponse(TrainingSet,1)=mod(sum(isingMachine,'all'),2);

    %Convert challenge to feature vector
    StateBitsPhi = Transform(neighbourMatrix, isingSize^2, maxconn); %
    challengebitsPhi(TrainingSet,:)= reshape(StateBitsPhi.',1,[]); %Convert to Linear 
end
toc

tic
fname = '.\';

filename= strcat(fname,'APUF_XOR_Challenge_Parity_64_1Million.csv');
writematrix([challengebitsPhi],filename)

filename= strcat(fname,'7-xorpuf_1M.csv');
writematrix(PUFResponse,filename)

fclose all;
toc