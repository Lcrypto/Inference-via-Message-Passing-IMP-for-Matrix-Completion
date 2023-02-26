function [Pr] = imp1(Re,Je,gu,gm)   
    nu = size(Re,1); nm = size(Re,2); 
    R2 = Re'; J2 = (R2>0); 
   
    % Generate VDVQ user class
    qvec = zeros(nm,1);
    for k=1:nm
        if(sum(Je(:,k))~=0)
           qvec(k) = sum(Re(:,k))/sum(Je(:,k));
        end
    end
    
    for i=1:log2(gu) 
        qvec = split_v01(qvec,2^i); class = ones(1,nm);
        for j=1:5
          [qvec,class] = vdlbgupdate_v01(R2,J2,qvec,@vq_meansq);
        end
    end
    classu = class; clear class
    
    % Generate VDVQ movei class
    qvec = zeros(nu,1);
    for k=1:nu
        if(sum(J2(:,k))~=0)
           qvec(k) = sum(R2(:,k))/sum(J2(:,k));
        end
    end
    
    for i=1:log2(gm) % modify index
        qvec = split_v01(qvec,2^i); class = ones(1,nu);
        for j=1:5
          [qvec,class] = vdlbgupdate_v01(Re,Je,qvec,@vq_meansq);
        end
    end
    classm = class; clear class

    % Convert to factor matrices
    eps = 1.5*10^(-1); 
    
    Ps = zeros(gu,size(Re,1));
    for i = 1:gu
        for j=1:nu
            if i == classu(j)
               Ps(i,j) = (1-eps);
            else
               Ps(i,j) = eps/(gu-1);
            end
        end
    end

    Pt = zeros(gm,size(Re,2));
    for i = 1:gm
        for j=1:nm
            if i == classm(j)
               Pt(i,j) = (1-eps);
            else
               Pt(i,j) = eps/(gm-1);
            end
        end
    end

    Pr = cell(1,5);
    for r=1:5
        Pr{r} = {zeros(gu,gm)};
    end
    for r=1:5
        for s=1:gu
            Ju = (classu==s);
            for t=1:gm
                Jm = (classm==t);
                Pr{r}{1}(s,t) = nnz(Re(Ju,Jm)==r); clear Jm
            end
            clear Ju
        end
    end

    dPr = zeros(gu,gm);
    for r=1:5
        dPr = dPr + Pr{r}{1};
    end
    dPr(dPr==0)=1;
    for r=1:5
        Pr{r}{1} = Pr{r}{1}./dPr;
    end
end


function [nqvec] = split_v01(qvec,n)
nvec1 = size(qvec,1);
nqvec = zeros(nvec1,n);

% split qvec into 2 times codebook
for i=1:n/2
    d = 0.2873*ones(nvec1,1);
    nqvec(:,i) = qvec(:,i) + d;
    nqvec(:,i+n/2) = qvec(:,i) - d;
end
end


function [qvec2,class,d] = vdlbgupdate_v01(train,train_flag,qvec,dist)
% First quantize to closest(class:codebook index, d:MS distortion)
[class,d] = vdvq(train,train_flag,qvec,dist); 
% Next, average over groups
qvec2 = zeros(size(qvec));
for i=1:size(qvec,2)                        
    D = sum(train(:,class==i),2);          
    N = sum(train_flag(:,class==i),2); 
    qvec2(N~=0,i) = D(N~=0)./N(N~=0);  
end
end


function [index,d] = vdvq(train,train_flag,qvec,dist)
% Quantize each to closest vector in qvec
index = zeros(1,size(train,2));              % initialization
d = zeros(1,size(train,2));                  % initialization

for i=1:size(train,2)                        % iterations #ntrain
  train3 = train(train_flag(:,i),i)';      % train vector without erasures
  qvec3 = qvec(train_flag(:,i),:)';          % codebook vectors corr.train2(optimize!)
  z = feval(dist,train3,qvec3);              % MS dist. between train(i,:)and codebook
  [m,mi] = min(z);                           % pick minimum among z   
  index(i) = mi;                             % index of the ith codebook
  d(i) = m;
end
end


function dist = vq_meansq(x,y)
% function dist = vq_meansq(x,y)
%
% Compute the mean squared distance between row vector x(=x(i,:))
% and set of row vectors y(i,:)
%
dist = sum(abs(y-x(ones(1,size(y,1)),:)).^2,2);  % pairwise sum along rows
end

