function []=digit2()

load digit.mat X T;

[dim, ~, numType] = size(X);
%numType = 3;

%�t�B�b�e�B���O
M = 5;
w = zeros(M, numType);
mu = zeros(dim, M, numType);
S = zeros(dim, dim, M, numType);
for i = 1:numType
    [w(:,i), mu(:,:,i), S(:,:,:,i)] = EM(M, X(:,:,i));
end

%�ΐ��ޓx�̌v�Z�i���Βl�j
[~,testcasesize,~] = size(T);
p = zeros(numType, testcasesize ,numType);
for expected = 1:numType
    t = T(:, :, expected);
    for i = 1:numType
        for j = 1:M
            muwi = mu(:, j,i);
            Swi = S(:,:,j,i);
            if cond(Swi) > 1e+15
                %���ٍs��ɂȂ�Ȃ��悤�ɔ��C��
                Swi = pinv(Swi, max(Swi,[],'all')*1e-10);
                p(expected,:,i) = p(expected,:,i)'...
                    + w(j,i)*(t'*Swi*muwi - muwi'*Swi*muwi/2);
            else
                p(expected,:,i) = p(expected,:,i)'...
                    + w(j,i)*(t'/Swi*muwi - muwi'/Swi*muwi/2);
            end
        end
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
disp(C);

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


end



