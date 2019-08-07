function [Xest] = admm_tv(Y, lmb, rho, dim,MAX_ITER)
    
    PRINT    = 1;
    TOL      = 1e-4;

    N        = dim(1);
    M        = dim(2);


    % Total Variational Operator
    Tmpn = zeros(N);
    Tmpm = zeros(M);    
    Tmpn(1,1) = -1;
    Tmpn(end,1) = 1;
    Tmpm(1,1) = -1;
    Tmpm(end,1) = 1;    
    Cn = real(ifft2(fft2(eye(N)).*fft2(Tmpn)));
    Cn(abs(Cn)<0.01) = 0;
    Cm = real(ifft2(fft2(eye(M)).*fft2(Tmpm)));
    Cm(abs(Cm)<0.01) = 0;    
    Lh = kron(speye(M), Cn);
    Lv = kron(Cm, speye(N));
 

    Lhinvl   = chol(Lh'*Lh + speye(N*M), 'lower');
    Lhinvu   = Lhinvl';
    Lvinvl   = chol(Lv'*Lv + speye(N*M), 'lower');
    Lvinvu   = Lvinvl';

    Z1       = zeros(N*M, 1);
    Z2       = zeros(N*M, 1);
    Z3       = zeros(N*M, 1);
    Z4       = zeros(N*M, 1);
    D1       = zeros(N*M, 1);
    D2       = zeros(N*M, 1);
    D3       = zeros(N*M, 1);
    D4       = zeros(N*M, 1);

    for itr = 1:MAX_ITER
        Z1_old = Z1;
        Z2_old = Z2;
        Z3_old = Z3;
        Z4_old = Z4;     
        
        % upadate f
        X      = (1/(1 + 2*rho))*( Y + rho*(Z1 - D1 + Z3 - D3) ); 
      
        % upadate x1
        Z1     = Lhinvu\(Lhinvl\(Lh'*(Z2 - D2) + (X + D1)));
 
        % upadate x2
        LhZ1   = Lh*Z1;
        Z2     = soft(LhZ1 + D2, lmb/rho);
        
        % upadate x3
        Z3     = Lvinvu\(Lvinvl\(Lv'*(Z4 - D4) + (X + D3)));
        
        % upadate x4
        LvZ3   = Lv*Z3;
        Z4     = soft(LvZ3 + D4, lmb/rho);
        
        % upadate d
        D1     = D1 + X - Z1;
        D2     = D2 + LhZ1 - Z2;
        D3     = D3 + X - Z3;
        D4     = D4 + LvZ3 - Z4;
        
        if mod(itr,10) == 1
            s = rho*norm(Z1-Z1_old + Z2-Z2_old + Z3-Z3_old + Z4-Z4_old,'fro');
            r = sqrt(norm(X-Z1,'fro')^2 + norm(LhZ1-Z2,'fro')^2 + norm(X-Z3,'fro')^2 ...
                                            + norm(LvZ3-Z4,'fro')^2);
                                       
            if r > 10*s
                rho = rho*2;
                D1 = D1/2;
                D2 = D2/2;
                D3 = D3/2;
                D4 = D4/2;
            elseif s > 10*r
                rho = rho/2;
                D1 = D1*2;
                D2 = D2*2;
                D3 = D3*2;
                D4 = D4*2;
            end                                        
            
            if  (s < TOL) && (r < TOL)
                break;
            end
            
            if PRINT  
%                 fprintf('itr = %f, res = %f, dual = %f, rho = %f\n', itr, r, s, rho)      
            end
        end
    end
    Xest = X;
end

function [x, res, iter] = cgsolve(A, b)

    tol = 1e-6; 
    maxiter = 20; 

    implicit = isa(A,'function_handle');

    x = zeros(length(b),1);
    r = b;
    d = r;
    delta = r'*r;
    delta0 = b'*b;
    numiter = 0;
    bestx = x;
    bestres = sqrt(delta/delta0); 
    while ((numiter < maxiter) && (delta > tol^2*delta0))

      % q = A*d
      if (implicit), q = A(d);  else  q = A*d;  end

      alpha = delta/(d'*q);
      x = x + alpha*d;

      if (mod(numiter+1,50) == 0)
        % r = b - Aux*x
        if (implicit), r = b - A(x);  else  r = b - A*x;  end
      else
        r = r - alpha*q;
      end

      deltaold = delta;
      delta = r'*r;
      beta = delta/deltaold;
      d = r + beta*d;
      numiter = numiter + 1;
      if (sqrt(delta/delta0) < bestres)
        bestx = x;
        bestres = sqrt(delta/delta0);
      end    

    end

    x = bestx;
    res = bestres;
    iter = numiter;
end

function z = soft(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end