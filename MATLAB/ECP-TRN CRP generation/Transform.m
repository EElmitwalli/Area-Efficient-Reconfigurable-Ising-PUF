function [APhi] = Transform(AChallenge, nRows, ChalSize )
% This code is modified from: https://github.com/scluconn/DA_PUF_Library
APhi = ones(nRows,ChalSize+1);
AChallenge = -2.*AChallenge+1;
    for j=1:ChalSize
        for k=j:ChalSize
            APhi(:,j) = APhi(:,j).*AChallenge(:,k);
        end
    end
end

