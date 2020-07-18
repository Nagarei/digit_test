function [w, mu, S] = EM(m, X)
%m個の正規分布でXをフィッティング

[dim, inputnum] = size(X);
%初期値
w = ones(m, 1)/m;
S = zeros(dim, dim, m);
Stmpmax = max(cov(X'),[],'all');
for i = 1:m
    vec = rand(dim)*Stmpmax;
    S(:,:,i) = vec*vec' + eye(dim);
end
mu = random('Normal', mean(X,'all'), abs(mean(X,'all'))/4, [dim m]);

L = -inf;
eta = zeros(inputnum, m);
while 1
    %E step
    for j = 1:m
        dx = X - repmat(mu(:,j), [1, inputnum]);
        Sj = S(:,:,j);
        eta(:,j) = w(j) * MyNormal(dx, Sj);
        %if cond(Sj) > 1e+15
        %    %特異行列にならないように微修正
        %    Sj = pinv(Sj, max(Sj,[],'all')*1e-10);
        %    eta(:,j) = w(j) * exp(-sum(dx'*Sj.*(dx'),2)/2);
        %else
        %    eta(:,j) = w(j) * exp(-sum(dx'/Sj.*(dx'),2)/2)...
        %        / (2*pi)^(dim/2) / sqrt(det(Sj));
        %end
    end
    etasum = sum(eta, 2);
    Lnew = sum(log(etasum));
    if Lnew-L<0.0001
        break;
    end
    L = Lnew;
    eta = eta ./ etasum;
    
    %M step
    wn = sum(eta,1);
    w = wn/inputnum;
    mu = (X*eta)./wn;
    for j=1:m
        dx = X - repmat(mu(:,j), [1, inputnum]);
        S(:,:,j) = dx*(dx'.*eta(:,j))/wn(j);
    end
end

end