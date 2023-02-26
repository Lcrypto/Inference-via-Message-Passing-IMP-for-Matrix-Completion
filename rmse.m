function [error] = rmse(Rt,Jt,Pr,Ps,Pt,gu,gm)
    Ghat = zeros(gu,gm);
    for r=1:5
        Ghat = Ghat + r*Pr{r}{1};
    end
    Rp = clip(Ps*Ghat*Pt').*Jt;
    error = full(sqrt(sum((Rt(:)-Rp(:)).^2)/nnz(Rt)));
end


function [R] = clip(R)
    R(R>5) = 5; R(R<1) = 1;
end

