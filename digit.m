function []=digit()

load digit.mat X T;

[dim, ~, numType] = size(X);
%numType = 3;

%平均値の計算
mu = zeros(dim, numType);
for i = 1:numType
    mu(:,i) = mean(X(:,:,i),2);
end

%共分散行列の計算
%カテゴリごとの共分散行列の算術平均
S = zeros(dim, dim);
S_rem = zeros(dim, dim);%誤差を減らすための一時変数
for i = 1:numType
    S(:,:) = float_add(S(:,:), S_rem(:,:), cov(X(:,:,i)'));
end
S(:,:) = S(:,:)/numType;

%対数尤度の計算（相対値）
[~,testcasesize,~] = size(T);
p = zeros(numType, testcasesize ,numType);
for expected = 1:numType
    t = T(:, :, expected);
    for i = 1:numType
        mui = mu(:,i);
        p(expected,:,i) = t'/S*mui - mui'/S*mui/2;
    end
end

%対数尤度基準で最大事後確率則を適用
[~, P] = max(p, [], 3);

%混同行列の生成
C = zeros(numType,numType);
for expected = 1:numType
    for i = 1:numType
        C(expected, i) = sum(P(expected,:)==i);
    end
end

%クラス別正誤率の出力
for i = 1:numType
    fprintf("class=%d: error=%f\n", rem(i,10), 1 - C(i,i) / sum(C(i,:)));
end

%recall,precisionの計算
recall = zeros(1, numType);
precision = zeros(1, numType);
F = zeros(1, numType);
for i = 1:numType
    precision(i) = C(i,i) / sum(C(:,i));
    recall(i) = C(i,i) / sum(C(i,:));
    F(i) = 2*recall(i)*precision(i)/(recall(i)+precision(i));
    fprintf("class=%d: recall=%f, precision=%f, F=%f\n",...
        rem(i,10), recall(i), precision(i), F(i));
end


%t = T(:, :, 2);

%invS = inv(S);
%p1 = t'/S*mu(:,1) - mu(:,1)'/S*mu(:,1)/2;
%p2 = t'/S*mu(:,2) - mu(:,2)'/S*mu(:,2)/2;

%result = sign(p1-p2);
%sum(result==1)
%sum(result==-1)

%t = T(:,find(result==-1), 1);
%t = T(:,69, 2);
%p1 = t'*invS*mu(:,1) - mu(:,1)'*invS*mu(:,1)/2;
%p2 = t'*invS*mu(:,2) - mu(:,2)'*invS*mu(:,2)/2;
%p1-p2
%colormap(gray);
%imagesc(reshape(t, [16 16])');

end
