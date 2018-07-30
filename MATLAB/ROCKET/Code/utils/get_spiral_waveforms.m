function [ktraj, G] = get_spiral_waveforms()

load ('Spirals_JonUoM_stack');
figure; plot(T,gx,'r'); hold on; plot(T,gy,'g'); plot(T,gz,'b');
xlabel('time (msec)'); ylabel('G/cm');