function [z,Relerrs] = PRSF(x,y, Params, A, At,z)    
u = Params.u0;
Relerrs =[];
for t = 1: Params.T,
    t
        grad = compute_grad(z, u, y, Params, A, At);
        z = z-Params.mu*grad;
        
%         ang = puma_ho(angle(z),.5,'verbose',verbose);
        
        Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro')]; % Relative error
        if norm(compute_grad(z,u,y,Params,A,At),'fro')< Params.y*u
            u = Params.y1*u;
        end
end