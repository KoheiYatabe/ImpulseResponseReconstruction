function z = RIRreconstADMM(data,lenX,obsIdx,preIdx,lambdaCoeff)
% ADMM algorithm for reconstructing impulse response from partial spectrum
%
% Only variables "data" and "obsIdx" are necessary for the ADMM algorithm.
% Other input variables are used for computing the weights "lambda" and "w"
% which can be constructed outside this function.
%
% Please refer to the following paper for details:
%
% Kohei Yatabe and Akiko Sugahara
% "Simulation of room impulse response using frequency-domain FEM and
%  convex-optimization-based post-processing"
% Applied Acoustics (2022)

maxIter = 3000;
mu = 1;
rho = 1.99;
[lambda,w] = weightGeneration(lenX,obsIdx,preIdx,lambdaCoeff);

X = zeros(lenX,1);
X(obsIdx) = data;
z = F(X);
u = 0;
for k = 1:maxIter
    x = proxPhi(z-u,mu*lambda);
    u = u + x - z;
    z = F(proxPsi(FH(x+u),data,obsIdx,mu*w));
    u = u + (rho-1)*(x - z);
end
z = z / sqrt(length(z));
end

function x = FH(x)
x = fft(x)/sqrt(length(x));
end

function x = F(x)
x = ifft(x,'symmetric')*sqrt(length(x));
end

function x = proxPhi(x,lambda)
x = max(1-lambda./abs(x),0) .* x;
end

function z = proxPsi(z,data,idx,w)
z = z./(1+2*w);
z(idx) = data;
end

function [lambda,w] = weightGeneration(lenX,obsIdx,preNum,lamdaCoeff)
lambda = lamdaCoeff*[ones(preNum,1); linspace(0,1,lenX-preNum)'.^2];
w = [linspace(1,0,min(obsIdx))'.^2; zeros(length(obsIdx)-2,1); linspace(0,1,floor(lenX/2)-max(obsIdx)+1)'.^2; ones(lenX-floor(lenX/2),1)];
end