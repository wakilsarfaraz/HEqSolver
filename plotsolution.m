%% Visualisation of the solution 
x = importdata('xcoordates.txt');
y = importdata('ycoordates.txt');
triang = importdata('triangles.txt');
U = importdata('SolutionMatrix.txt');
T = importdata('timesteps.txt');


v = VideoWriter ('DiffusionEquation.avi');
open (v) ;

for j = 1 : length(U(:,1))
    trisurf(triang, x, y, U(j,:))
    colorbar 
    colormap(hsv)
    grid on
    shading interp
    zlim([0 1])
    title(['Diffusion of heat at time t = ' num2str(T(j))],'fontsize',20)
    xlabel('x','fontsize',16)
    ylabel('y','fontsize',16)
    zlabel(['U(x,y,',num2str(T(j)),')'],'fontsize',16)
    frame = getframe (gcf);
    writeVideo (v, frame); 
    pause(1e-5)
end

close(v)

  set(findobj('type','legend'),'fontsize',18)
set(findobj('type','axes'),'fontsize',18)
