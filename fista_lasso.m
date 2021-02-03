function x=fista_lasso(A,b,mu,maxiter)
    tic
    [m,n]=size(A);
    ATA=A'*A;
    L=1*eigs(ATA,1);
    %% Calculate Prox_{\gamma g}(y)
    function value=Prox(y,mu,gamma)
        value=[];
        for i=1:length(y)
            if y(i)>mu*gamma
                value=[value;y(i)-mu*gamma];
            elseif y(i)<-mu*gamma
                value=[value;y(i)+mu*gamma];
            else 
                value=[value;0];
            end
        end
    end   
    %% Calculate PL function
    x=zeros(n,1);
    y=x;
    function valuepl=PL(y)
        valuepl=Prox(y-(1/L)*A'*(A*y-b),mu,1/L);
    end
    steps=0;
    t=1;
    while steps<maxiter
        x0=x;
        x=PL(y);
        t0=t;
        t=0.5*(1+sqrt(1+4*(t^2)));
        y=x+((t0-1)/(t))*(x-x0);
        steps=steps+1;
        if rem(steps,100)==0
            fprintf('iter: %d\n',steps)
            fprintf('Value_Fista: %10.14f\n',0.5*norm(A*x-b)^2+mu*norm(x,1))
            fprintf('time: %d\n',toc)
        end
    end
end