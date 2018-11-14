function ktraj_plotter(ktrajs)

for proj =1:size(ktrajs,1)
    
    if(size(ktrajs,3) == 3)
        figure(1003);
        plot3(squeeze(ktrajs(proj,:,1)), squeeze(ktrajs(proj,:,2)),squeeze(ktrajs(proj,:,3)), 'r*');
        hold on; xlabel('kx - /mm'); ylabel('ky- /mm'); zlabel('kz - /mm'); grid on;
    else
    end
    
    
    
end
