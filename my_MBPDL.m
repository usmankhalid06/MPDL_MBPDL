function [D,X,Err,C]= my_Proposed(Y,Dp,Xp,K,spa,l1,l2,l3,l4,nIter,fac,TC,SM)
    A = eye(size(Dp,2),K);
%     D = Dp*A;  
    D = Y(:,1:K);   D = D*diag(1./sqrt(sum(D.*D))); 
    B = zeros(K,size(Xp,1));
    X = zeros(size(D,2),size(Y,2)); 
    fprintf('Iteration:     ');
    for iter=1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
        Dold = D;

        X2 = abs(X);
        for jjj=1:size(X2,1)
            X2(jjj,:) =(X2(jjj,:) - min(X2(jjj,:))) / ( max(X2(jjj,:)) - min(X2(jjj,:)) );
        end
        for jjj=1:fac
            [~,ind4(jjj)] = max(abs(corr(Xp(jjj,:)',X2')));
        end


        for j =1:size(D,2)
            X(j,:) = 0; A(:,j) = 0; B(j,:) = 0;
            E = Y-D*X;


            if any(j==ind4)
                lam1 = l1;  lam2= l2;  
            else
                lam1 = l3;  lam2 = l4;  
            end

            xk = D(:,j)'*E;
            thr = spa./abs(xk);
            xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
            [~,bb1]= sort(abs(Xp*xkk'),'descend');
            ind1 = bb1(1:lam2);
            B(j,ind1)= xkk*pinv(Xp(ind1,:));  %xkk*Xp(ind1,:)'/(Xp(ind1,:)*Xp(ind1,:)'); %
%             B(j,:) = B(j,:) / sum(B(j,:));
            X(j,:) = B(j,:)*Xp;

            rInd = find(X(j,:));
            if (length(rInd)<1)

            else
                tmp3 = E(:,rInd)*X(j,rInd)';
                [~,bb]= sort(abs(Dp'*tmp3),'descend');
                ind = bb(1:lam1);
                A(ind,j)= (Dp(:,ind)'*Dp(:,ind))\Dp(:,ind)'*tmp3;
%                 A(:,j) = A(:,j) / sum(A(:,j));
                A(:,j) = A(:,j)./norm(Dp*A(:,j));
                D(:,j) = Dp*A(:,j);
            end

        end      
        Err(iter) = sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold));
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
end



