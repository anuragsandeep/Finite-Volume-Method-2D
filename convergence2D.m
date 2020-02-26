%% Heat Transfer in 2-D using line by line visitation
%% with Gauss-Seidel and TDMA 
%% Anurag Sandeep K. 
clc
clear all

%% CHECK FOR CONVERGENCE
alpha=1.0; %1.0:0.03:1.39; % over-relaxation factor
NX=[17 27 37 51]; % ITMAX+2 in x
NY=NX;            % JTMAX+2 in y
ITMX=NX-2;
TCentre=[];
 for ia=1:length(NX)
     [T{ia},ITER(ia)]=heatTransfer2D(NX(ia),NY(ia),alpha);
 end

 for ia=1:length(NX)
     TCenter(ia) = T{ia}((NX(ia)+1)/2,(NY(ia)+1)/2);
 end
 
  
% P L O T T I N G

% plot(alpha,ITER,'LineWidth',2); hold on
% plot(1.36,87,'MarkerFaceColor',[1 0 0],'Marker','o','Color',[1 0 0])
% xlabel('$\omega (over-relaxation factor)$','Interpreter','Latex')
% xlabel('$\omega$ (over-relaxation factor)','Interpreter','Latex')
% ylabel('Iterations')
% title('Iterations for convergence')
% text(1.3,400,'\omega_{opt}=1.36')

% % For contour plot
% [M,N] = size(T{1});
% [x,y] = meshgrid(1:N,1:M); 
% surf(x,y,T{1}'); shading interp
% view(2); colormap('hot');colorbar
