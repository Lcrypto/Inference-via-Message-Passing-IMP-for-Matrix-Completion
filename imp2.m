function [Pt Ps] = imp2(Re,Pr,gu,gm,max_iter)

    GMatrix = zeros(gu,gm,5);
    for i=1:5
        GMatrix(:,:,i) = Pr{i}{1};
    end

    % Generate graph
    [MOV USER] = generate_graph(Re);
    MOV(MOV==0) = -1;
    USER(USER==0) = -1;

    % Get max user and movie degrees
    max_movie_degree = size(MOV,2);
    max_user_degree = size(USER,2);

    % Initialize messages as uniform
    nu = size(Re,1); nm = size(Re,2); 
    mnode_init = 1/gm*ones(nm,gm);
    unode_init = 1/gu*ones(nu,gu);
    % movie user gest
    [Pt Ps] = mpa_bk(max_iter,gu,MOV-1,max_movie_degree,USER-1, ...
              max_user_degree,mnode_init,unode_init, ...
              GMatrix(:),Re,0);
end

