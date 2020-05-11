function S = Entropy(rho)
ew=eig(rho);
ew=ew.*heaviside(ew);ew(ew==0)=[];
S=-real(sum(ew.*log(ew)));
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

