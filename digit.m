function []=digit()

load digit.mat X T;

[dim, ~, numType] = size(X);
%numType = 3;

%���ϒl�̌v�Z
mu = zeros(dim, numType);
for i = 1:numType
    mu(:,i) = mean(X(:,:,i),2);
end

%�����U�s��̌v�Z
%�J�e�S�����Ƃ̋����U�s��̎Z�p����
S = zeros(dim, dim);
S_rem = zeros(dim, dim);%�덷�����炷���߂̈ꎞ�ϐ�
for i = 1:numType
    S(:,:) = float_add(S(:,:), S_rem(:,:), cov(X(:,:,i)'));
end
S(:,:) = S(:,:)/numType;

%�ΐ��ޓx�̌v�Z�i���Βl�j
[~,testcasesize,~] = size(T);
p = zeros(numType, testcasesize ,numType);
for expected = 1:numType
    t = T(:, :, expected);
    for i = 1:numType
        mui = mu(:,i);
        p(expected,:,i) = t'/S*mui - mui'/S*mui/2;
    end
end

%�ΐ��ޓx��ōő厖��m������K�p
[~, P] = max(p, [], 3);

%�����s��̐���
C = zeros(numType,numType);
for expected = 1:numType
    for i = 1:numType
        C(expected, i) = sum(P(expected,:)==i);
    end
end

%�N���X�ʐ��뗦�̏o��
for i = 1:numType
    fprintf("class=%d: error=%f\n", rem(i,10), 1 - C(i,i) / sum(C(i,:)));
end

%recall,precision�̌v�Z
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
