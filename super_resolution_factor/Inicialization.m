function [z,Relerrs] = Inicialization(x,y, Params, A, At)    
%% Initialization
    const = Params.const;
    npower_iter = Params.npower_iter;                                      % Number of power iterations
    
    z0 = randn(Params.n1,Params.n2); 
    z0 = z0/norm(z0,'fro');                                                % Initial guess
    itervalo=Params.itervalo;
    normest = sqrt(sum(y(:).^2)/numel(y(:))/const^2);                                 % Estimate norm to scale eigenvector
        
    [~,B] = sort(y(:).^2,'descend');
    ytr = zeros(size(y));
    ytr(B(1:ceil(Params.m/Params.p))) = 1;
    
   filt = zeros(Params.n1,Params.n2);
   filt(itervalo,itervalo) =1;
   se = fspecial('gaussian',5);
    for tt = 1:npower_iter,                                                % Truncated power iterations
        tt
        z0 = At(ytr.* A(z0))./(Params.m*ceil(Params.m/Params.p));
        z0 = filt.*z0;
%         z0 = imfilter(z0,se);
        z0 = z0/norm(z0,'fro');
%         ang = puma_ho(angle(z0(193:320,193:320)),.5,'verbose','no');
%         z(193:320,193:320) = abs(z0(193:320,193:320)).*exp(1i*ang);
    end
       
    z = normest * z0;                   % Apply scaling
    u = Params.u0;
    Relerrs = norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro'); % Initial rel. error    