function phi= evolution(I,phi0,mu,nu,lambda,delta_t,epsilon,numIter)
%   new_revolution(I,phi0,g,gx,gy,mu,nu,lambda,delta_t,epsilon,numIter)
%   input: 
%       I: input image
%       phi0: level set function to be updated
%       g: edge detector g
%       gx,gy: gradient of edge detector g
%       mu: distRictTerm coefficient
%       nu: weight for area term, default value 0
%       lambda: lengthTerm coefficient
%       delta_t: time step
%       epsilon: parameter for computing smooth Heaviside and dirac function
%       numIter: number of iterations
%   output: 
%       phi: updated level set function
sigma = 1.5;             %   选取论文给定参数
g=edge_detector(I,sigma);
[gx,gy]=gradient(g);
% I = BoundMirrorExpand(I); % 镜像边缘延拓
phi = BoundMirrorExpand(phi0);

% gaussian filtering
[px,py]=gradient(phi);
pSum=sqrt(px.^2 + py.^2 + 1e-10);
px=px./pSum;
py=py./pSum;

for k = 1 : numIter
    phi = BoundMirrorEnsure(phi);
    delta_h = Delta(phi,epsilon);
    Curv = curvature(phi);
    phi = phi + delta_t * (mu * (4 * del2(phi) - Curv) + lambda * delta_h .*( px.*gx + py.*gy + g.*Curv) + nu * g .* delta_h);
end
phi = BoundMirrorShrink(phi);



