RMS(X,Y,Z) = sqrt((vec(X).-mean(vec(X)))^2+(vec(Y).-mean(vec(Y)))^2+(vec(Z).-mean(vec(Z)))^2)
RMS(X,Y) = sqrt((vec(X).-mean(vec(X)))^2+(vec(Y).-mean(vec(Y)))^2)
RMS(X) = sqrt(vec(X)-mean(vec(X).^2))