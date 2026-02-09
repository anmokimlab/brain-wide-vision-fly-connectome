function [xstar, ystar, zstar, info] = closest_point_poly22(P, p, opts)
    if nargin<3, opts=struct; end
    if ~isfield(opts,'tol'),   opts.tol = 1e-16; end
    if ~isfield(opts,'maxit'), opts.maxit = 100;   end
    p00=P(1); p10=P(2); p01=P(3); p20=P(4); p11=P(5); p02=P(6);
    px=p(1); py=p(2); pz=p(3);
    f  = @(x,y) p00 + p10*x + p01*y + p20*x.^2 + p11*x.*y + p02*y.^2;
    fx = @(x,y) p10 + 2*p20*x + p11*y;
    fy = @(x,y) p01 + p11*x + 2*p02*y;
    fxx=2*p20; fxy=p11; fyy=2*p02;

    x=px; y=py;   % init: vertical drop
    for k=1:opts.maxit
        F=f(x,y); gx=fx(x,y); gy=fy(x,y); dz=F-pz;
        r1 = x - px + dz*gx;  r2 = y - py + dz*gy;
        if hypot(r1,r2) < opts.tol*(1+hypot(px,py)), break; end
        J11 = 1 + gx*gx + dz*fxx;   J12 = gx*gy + dz*fxy;
        J22 = 1 + gy*gy + dz*fyy;   J21 = J12;
        delta = -[J11 J12; J21 J22] \ [r1; r2];

        % 간단 백트래킹
        step=1; old=hypot(r1,r2);
        for bt=1:6
            xn=x+step*delta(1); yn=y+step*delta(2);
            Fn=f(xn,yn); gxn=fx(xn,yn); gyn=fy(xn,yn); dzn=Fn-pz;
            new=hypot(xn-px + dzn*gxn, yn-py + dzn*gyn);
            if new <= 0.9*old || step < 1/64, x=xn; y=yn; break; end
            step=step*0.5;
        end
    end
    xstar=x; ystar=y; zstar=f(x,y);
    info.iter=k; info.dist=norm([x-px, y-py, f(x,y)-pz]);
end
