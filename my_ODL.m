function [D, X, Err,C] = my_ODL(YY, L, nIter,lambda,bs,nIter2,TC,SM)
%     D = dctbases(size(YY,1),L);
    D = YY(:,1:L);D = D*diag(1./sqrt(sum(D.*D))); 
    param.mode = 2; 
    param.lambda = lambda; 
    F = 0; 
    E = 0;
    fprintf('Iteration:     ');
    for iter =1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
%         Y = YY(:,bs*(iter-1)+1:bs*iter);
        Y = YY; 
        Dold = D;
%         X = mexLasso(Y, D, param);
        X = customLasso(Y, D, lambda, nIter);
        F = F+X*X';
        E = E+Y*X';
        Dp = D;
        iter2 = 0;
        while (iter2 < nIter2)
            iter2 = iter2 + 1;
            for k = 1: size(D,2)
                if(F(k,k) ~= 0)
                    tmpD = 1.0/F(k,k) * (E(:,k) - D*F(:, k)) + D(:,k);
                    D(:,k) = tmpD/(max( norm(tmpD,2),1));
                end
            end
            if (norm(D - Dp, 'fro')/numel(D) <1e-9)
                break;
            end
            Dp = D;
        end
        Err(iter)= sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold));
        K2 = size(TC,2);
        [~,~,ind]=sort_TSandSM_spatial(TC,SM,D,X,K2);
        for ii =1:K2
            TCcorr(ii) =abs(corr(TC(:,ii),D(:,ind(ii))));
            SMcorr(ii) =abs(corr(SM(ii,:)',X(ind(ii),:)'));
        end
        cTC = sum(TCcorr');
        cSM = sum(SMcorr');
        C(iter) =cTC+cSM;             
    end
%     X = mexLasso(YY, D, param);
    X = customLasso(Y, D, lambda, nIter);

end


function X = customLasso(Y, D, lambda, numIter)
    [~, n] = size(D);
    [~, p] = size(Y);
    X = zeros(n, p);
    
    for j = 1:p
        x = zeros(n, 1); % Initial solution for each column of Y
        for iter = 1:numIter
            for i = 1:n
                % Compute the residual without the contribution of the i-th basis
                r = Y(:, j) - D * x + D(:, i) * x(i);
                
                % Soft thresholding
                a = D(:, i)' * r;
                x(i) = sign(a) * max(abs(a) - lambda, 0) / (D(:, i)' * D(:, i));
            end
        end
        X(:, j) = x;
    end
end


% function [D, X, Err] = my_ODL(Y, K, lambda, nIter)
%     Di = dctbases(size(Y,1),K);
%     X = zeros(K size(Y,2));
% 	iter = 0;
%     param.mode = 2; 
%     param.lambda = lambda; 
%     fprintf('Iteration:     ');
%     sizeD = numel(D);    
% 	while (iter < nIter)
%         Dold = D;
% 		iter = iter + 1;
%         fprintf('\b\b\b\b\b%5i',iter);
%         X = mexLasso(Y, D, param);
% 		
%         F = X*X'; E = Y*X';
%         Dp = D;
%         iiter = 0;
%         while (iiter < nIter)
%             iiter = iiter + 1;
%             for i = 1: size(D,2)
%                 if(F(i,i) ~= 0)
%                     a = 1.0/F(i,i) * (E(:,i) - D*F(:, i)) + D(:,i);
%                     D(:,i) = a/(max( norm(a,2),1));
%                 end
%             end
%             %% check stop condition
%             if (norm(D - Dp, 'fro')/sizeD < 1e-6)
%                 break;
%             end
%             Dp = D;
%         end
%         Err(iter)= sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold));
%     end
%     
% end
