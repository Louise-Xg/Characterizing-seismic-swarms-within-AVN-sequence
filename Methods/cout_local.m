@author: David Marsan

function [J,Mmodele, Mmodele_liss] = cout_local(theta,t,M,m,I,dn)
% modèle : M = -a * log(t-t_MS) + m_MS - b, et M = c si M<c. Si plusieurs
% MS : consièdre celui donnant le max de M.

N=length(I);
a=theta(1:N); b=theta(N+1:2*N); c=theta(end);
Mmodele=M*0+c; 
for i=1:N
    n=I(i);
    J=find(t>t(n));
    tmp=-a(i)*log(t(J)-t(n))+m(n)-b(i);
    Mmodele(J)=max([Mmodele(J) ; tmp]);
end
% -> lissage de Mmodele
Mmodele_liss = movmean(Mmodele,dn);
J=mean((M-Mmodele_liss).^2); %%rps l'écart entre le modèle et les données observées <=> erreur

end
