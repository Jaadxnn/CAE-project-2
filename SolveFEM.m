function [U, Ux, Uy, Umag, ...
          Ex, Ey, Gxy, ...
          Sx, Sy, Sxy] = ...
          SolveFEM( ...
          n_nodes, n_elements, ...
          ncon, X, Y, ...
          E, v, t, ...
          F, dzero, NDU)



            K = zeros(2*n_nodes);
            for i = 1:n_elements

                n1 = ncon(i,1);
                n2 = ncon(i,2);
                n3 = ncon(i,3);
            
                x1 = X(n1); y1 = Y(n1);
                x2 = X(n2); y2 = Y(n2);
                x3 = X(n3); y3 = Y(n3);
            
                Ae = 0.5 * det([1 x1 y1;
                                1 x2 y2;
                                1 x3 y3]);
            
                Ae = abs(Ae);
            
                b1 = y2 - y3;
                b2 = y3 - y1;
                b3 = y1 - y2;
            
                c1 = x3 - x2;
                c2 = x1 - x3;
                c3 = x2 - x1;
            
                B = (1/(2*Ae))*[b1 0  b2 0  b3 0;
                                 0 c1 0 c2 0 c3;
                                 c1 b1 c2 b2 c3 b3];
            
                D = (E/(1-v^2))*[1 v 0;
                                 v 1 0;
                                 0 0 (1-v)/2];
%Checking
disp(size(B))
disp(size(D))
disp(E)
disp(v)
disp(t)
                KE = (t)*Ae*(B.')*D*B;
            
                ROC = [2*n1-1 2*n1 2*n2-1 2*n2 2*n3-1 2*n3];

                for ix = 1:6
                    for jx = 1:6
                        K(ROC(ix),ROC(jx)) = K(ROC(ix),ROC(jx)) + KE(ix,jx);
                    end
                end
            end

%% =========================
% 8. APPLY BCs
% =========================

KM = K;

for k = 1:NDU
    n = dzero(k);
    KM(n,:) = 0;
    KM(:,n) = 0;
    KM(n,n) = 1;
end

%% =========================
% 9. SOLVE
% =========================

U = KM \ F;
U = U;

%*********POST PROCCESSING**************
Ex = zeros(n_elements,1);
Ey = zeros(n_elements,1);
Gxy = zeros(n_elements,1);

Sx = zeros(n_elements,1);
Sy = zeros(n_elements,1);
Sxy = zeros(n_elements,1);

Ae_all = zeros(n_elements,1);

for i = 1:n_elements

    n1 = ncon(i,1);
    n2 = ncon(i,2);
    n3 = ncon(i,3);

    x1 = X(n1); y1 = Y(n1);
    x2 = X(n2); y2 = Y(n2);
    x3 = X(n3); y3 = Y(n3);

    Ae = 0.5 * det([1 x1 y1;
                    1 x2 y2;
                    1 x3 y3]);

    Ae = abs(Ae);
    Ae_all(i) = Ae;

    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;

    B = (1/(2*Ae))*[b1 0  b2 0  b3 0;
                     0 c1 0 c2 0 c3;
                     c1 b1 c2 b2 c3 b3];

    D = (E/(1-v^2))*[1 v 0;
                     v 1 0;
                     0 0 (1-v)/2];

    d = [U(2*n1-1); U(2*n1);
         U(2*n2-1); U(2*n2);
         U(2*n3-1); U(2*n3)];

    e = B*d;
    Sigma = D * e;

    Ex(i) = e(1);
    Ey(i) = e(2);
    Gxy(i) = e(3);

    Sx(i)  = Sigma(1);
    Sy(i)  = Sigma(2);
    Sxy(i) = Sigma(3);
end

Ex = Ex;
Ey = Ey;
Gxy = Gxy;

Sx = Sx;
Sy = Sy;
Sxy = Sxy;

Ux = U(1:2:end);
Uy = U(2:2:end);
Umag = sqrt(Ux.^2 + Uy.^2);

end