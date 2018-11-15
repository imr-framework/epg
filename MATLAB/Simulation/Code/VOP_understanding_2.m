%% 3D visualization
close all;
for k=1:72
    ind_int = find(VOP_ind ==k);
    if(~isempty(ind_int))
%     figure;
    plot(cv_field(ind_int),xshim_LSAR(ind_int),'r*');axis([1 3 0 15]);
    pause;
    end
end
%%
figure;
plot3(VOP_ind,abs(xshim_LSAR),mls_error,'r*');
grid on;xlabel('VOP index');ylabel('Local SAR'); zlabel('MLS error');


%% 
figure;
subplot(131);plot(VOP_ind,abs(xshim_LSAR),'r*');xlabel('VOP index');ylabel('L-SAR');
subplot(132);plot(VOP_ind,abs(mls_error),'r*');xlabel('VOP index');ylabel('MLSE');
subplot(133);plot(mls_error,abs(xshim_LSAR),'r*');xlabel('MLSE');ylabel('L-SAR');


%%
figure;
subplot(121);plot(mls_error,abs(xshim_LSAR),'r*');xlabel('MLSE');ylabel('L-SAR');
subplot(122);plot(abs(cv_field),abs(xshim_LSAR),'r*');xlabel('CV field');ylabel('L-SAR');

%%
figure;plot(VOP_ind,abs(cv_field),'r*');xlabel('VOP index');ylabel('CV');
plot(VOP_ind,abs(mls_error),'r*');xlabel('VOP index');ylabel('MLSE');axis([1 72 0 3.5])

%% Filtered solutions
ind_mlse = find(mls_error < 3.5);
% ind_mlse = find(xshim_LSAR < 2);

%%
figure;
plot(ind_mlse,mls_error(ind_mlse),'r*-');

figure;
plot3(VOP_ind(ind_mlse),abs(xshim_LSAR(ind_mlse)),mls_error(ind_mlse),'r*');
grid on;xlabel('VOP index');ylabel('Local SAR'); zlabel('MLS error');

%%
figure;
subplot(131);plot(VOP_ind(ind_mlse),abs(xshim_LSAR(ind_mlse)),'r*');xlabel('VOP index');ylabel('L-SAR');
subplot(132);plot(VOP_ind(ind_mlse),abs(mls_error(ind_mlse)),'r*');xlabel('VOP index');ylabel('MLSE');
subplot(133);plot(mls_error(ind_mlse),abs(xshim_LSAR(ind_mlse)),'r*');xlabel('MLSE');ylabel('L-SAR');

%% 


