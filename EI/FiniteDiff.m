clear
close all

nx = 50; ny = 50;

G = sparse(nx*ny, nx*ny);

a = 0;

num = 16;

for i = 0:ny-1
    if i == 0 || i == ny - 1
        for j = 0:nx-1
            G(ny*i+j+1,ny*i+j+1) = 1;
        end
    else
        G(i*ny+1, i*ny+1) = 1;
        G(i*ny+nx, i*ny+nx) = 1;
    end
end

for i = 0:ny
    if ((0 < i) && (i < (ny - 1)))
        for j = 0:nx
            if ((0 < j) && (j < (nx - 1)))
                G(i*ny+j + 1, i*ny+j + 1) = -4 - a; %Center
                G(i*ny+j + 1, i*ny+j+1 + 1) = 1;     %Right
                G(i*ny+j + 1, i*ny+j-1 + 1) = 1;     %Left
                G(i*ny+j + 1, i*ny+j+nx + 1) = 1;    %Below
                G(i*ny+j + 1, i*ny+j-nx + 1) = 1 ;   %Above
            end
        end
    end
end

spy(G);

[E,D] = eigs(G,num,'SM')

spy(D);


for k = 1:num
    subplot(4,4,k);
    surf(reshape(E(:,k),nx,ny))
end