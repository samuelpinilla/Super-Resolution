function grad = compute_grad_super_v2(z, u, y, Params, A, At,x1)
    m = Params.m;
    x1 = reshape(x1,size(z));
    yz = sqrt(abs(A(z)).^2+u^2);
    grad = 2/m*At(((yz-y)./yz).*A(z))+1.*(z-x1);
%     norm(grad,'fro')
end

  