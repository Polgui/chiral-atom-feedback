function rhovec=mat2vec(rho)
dims=size(rho);
rhovec=reshape(rho,[dims(1)*dims(2),1]);
end