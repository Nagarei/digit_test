function result = MyNormal(dx, sigma)

[dim, testcases] = size(dx);
if cond(sigma) < 1e+15 && det(sigma) > 1e-10
    result =  exp(-sum(dx'/sigma.*(dx'),2)/2)...
        / (2*pi)^(dim/2) / sqrt(det(sigma));
    return;
end

% fprintf("MYNORMAL %f %f", cond(sigma), det(sigma))

%ì¡àŸçsóÒÇ…Ç»ÇÁÇ»Ç¢ÇÊÇ§Ç…èCê≥
[P, L] = eig(sigma);
dxP = P'*dx;
okmin = max(L,[],'all')*1e-6;
zerocount = 0;
for d = 1:dim
    if L(d,d) < okmin
        zerocount = zerocount + 1;
    end
end
convert = zeros(dim, dim - zerocount); %convert*LÇ≈ècÇÃ0ê¨ï™ÇèúÇ≠
convert_size = 0;
ng_arr = ones(1, testcases);
for d = 1:dim
    if L(d,d) < okmin
        ng_arr = ng_arr .* (abs(dxP(d,:)) <= 10*L(d,d));
    else
        convert_size = convert_size+1;
        convert(d, convert_size) = 1;
    end
end

dx2 = convert'*dxP;
sigma2 = convert'*L*convert;
result =  exp(-sum(dx2'/sigma2.*(dx2'),2)/2)...
    / (det((2*pi)^(dim/2) * sqrt(sigma2)));
%result = result .* ng_arr';

return;
    
end
