function rho=vec2mat(rhovec)
dimvec=size(rhovec);
dim=sqrt(dimvec(1));
rho=reshape(rhovec,dim,dim);
end