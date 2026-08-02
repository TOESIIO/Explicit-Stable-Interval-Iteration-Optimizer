%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE (MODIFIED FROM top88) %%%%
function top88ESIIO_H_ESIIO(nelx,nely,volfrac,penal,rmin,ft)
nelx=800; nely=500; volfrac=0.4; penal=3;  rmin=3; ft=1;
%% MATERIAL PROPERTIES
k0 = 1;
kmin = 1e-3;
%% PREPARE FINITE ELEMENT ANALYSIS
KE=[2/3 -1/6 -1/3 -1/6
    -1/6 2/3 -1/6 -1/3
    -1/3 -1/6 2/3 -1/6
    -1/6 -1/3 -1/6 2/3];
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
edofMat = repmat(edofVec,1,4)+repmat([0 1 nely+[1 2]],nelx*nely,1);
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);
U = zeros((nely+1)*(nelx+1),1);
fixeddofs = [1:(nely+1)];
alldofs = [1:(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
F = sparse(1:(nely+1)*(nelx+1),1,1/((nely+1)*(nelx+1)),(nely+1)*(nelx+1),1);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
loop = 0;
change = 1;
p = 10; lam2 = 1;
%% START ITERATION
while loop < 200 && change > 0.01
    loop = loop + 1;
    xPhys = reshape((ft==1)*x(:) + (ft==2)*(H*x(:))./Hs, nely, nelx);
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(kmin+xPhys(:)'.^penal*(k0-kmin)),16*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((kmin+xPhys.^penal*(k0-kmin)).*ce));
    dc = -penal*(k0-kmin)*xPhys.^(penal-1).*ce;
    v = mean(xPhys(:)) - volfrac;
    dv = ones(nely,nelx) / (nelx*nely);
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    move = 0.2;  m = 1; n = nelx*nely; v = mean(xPhys(:))-volfrac;
    xmax = min(1,x(:)+move); xmin = max(1e-6,x(:)-move);
    if loop == 1
        lam1=zeros(1,m); mu=10*ones(1,m); fold=v; f1=0; f2=0; cc=0;
    end
    [xnew, lam1, mu, lam2, fold, f1, f2, cc] = ESIIO(...
            m, n, loop, x(:), xmin, xmax, c, dc(:), v, dv(:)', move, p,...
            lam1, mu, lam2, fold, f1, f2, cc);
    if ~mod(loop, 20) && loop >= 40, p = min(20*p, 2e4); end
    %% CONVERGENCE JUDGMENT
    change = max(abs(xnew(:)-x(:)));
    x = reshape(xnew,nely,nelx);
    %% PRINT RESULTS
    fprintf(' It.:%5i Comp:%11.4f Vol.:%7.3f ch.:%7.3f\n',...
        loop, c, mean(xPhys(:)), change);
    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end