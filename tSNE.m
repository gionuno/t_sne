function [Y,S,P,E] = tSNE(X,D,H,T)
    N = size(X,1);
    M = size(X,2);
    
    E = zeros(T,1);
    
    Y = 1e-4*randn(N,D);
    mY = zeros(N,D);
    nY = zeros(N,D);
    
    mu = 0.5;
    nu = 0.5;
    
    S = ones(N,1);
    
    for i = 1:N
        disp(i);
        ls = -10.0;
        hs =  10.0;
        pl = zeros(N,1);
        for j = 1:N
            pl(j) = exp(-0.5*((2^ls)*norm(X(i,:)-X(j,:))/M)^2);
        end
        pl(i) = 1e-100;
        pl = pl / sum(pl);
        el = -pl'*log2(pl)-log2(H);
        while abs(ls-hs)>1e-4            
            ms = 0.5*(ls+hs);
            pm = zeros(N,1);
            for j = 1:N
            	pm(j) = exp(-0.5*((2^ms)*norm(X(i,:)-X(j,:))/M)^2);
            end
            pm(i) = 1e-100;
            pm = pm / sum(pm);
            em = -pm'*log2(pm)-log2(H);
            if abs(em) < 1e-4
                ls = ms;
                hs = ms;
                break;
            end
            if sign(em)*sign(el) > 0 
                ls = ms;
                el = em;
            else
                hs = ms;                
            end
        end
        S(i) = 2^(0.5*(hs+ls))/M;
        disp(S(i));
    end
    
    P = zeros(N,N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                P(i,j) = exp(-0.5*(S(i)*norm(X(i,:)-X(j,:)))^2);
            end
        end
        P(i,:) = P(i,:)/sum(P(i,:));
    end
    P = 0.5*(P+P')/N;
    
    C = zeros(N,N);
    for t = 1:T
        disp(t);
        for i = 1:N
            for j = (i+1):N
                C(i,j) = 1.0/(1.0+norm(Y(i,:)-Y(j,:))^2);
                C(j,i) = C(i,j);
            end
        end
        Q = C / sum(sum(C));
        E(t) = -mean(sum(P.*log(Q+1e-10)));
        disp(E(t));
        A = (P-Q).*C;
        M = diag(sum(A,2))-A;
        dY = M*Y;
        mY = mu*mY + (1-mu)*dY;        
        nY = nu*nY + (1-nu)*dY.^2;

        mY_ = mY / (1-mu^t);
        nY_ = nY / (1-nu^t);
        
        Y = Y - 5e-1* mY_ ./ (sqrt(nY_)+1e-8);
    end
end