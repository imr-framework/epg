%------------------------------------------------------------------------
% df_read:  18-10-07. Reads Philips DF files
% SJM.      USAGE: image = df_read('filename');
%------------------------------------------------------------------------

function im = df_read(filename)

% 22-11-10: only add .df if not present
if isempty(strfind(filename,'.df'))
    fname = [filename '.df'];
else
    fname=filename;
end
id = fopen(fname,'r','l');

% get header info
hd = fread(id,4,'uint32');

if hd(1) ~= hex2dec('DEADFACE')
    disp('error, not df format');
    return
end
    
% hd(2)/hd(3) = FE/PE numbers, hd(4) = num slices. factor 2 = real/imag
im = fread(id,prod(hd(2:4))*2,'*float');

% reshape. looks like the complex numbers are split up so that it is real
% then imag. but i'm not sure about the order of slice/cplx

if hd(4)==1
    if(numel(im) ~= prod(hd(2:3)))
        im = reshape(im,[hd(2) hd(3) 2]);
        im = complex(im(:,:,1),im(:,:,2));
    else
        im = reshape(im,[hd(2) hd(3)]);
    end
else
    %% added code to error handle BodySegmentation.df for 1 channel data
if(numel(im) == prod(hd(2:4)))
        im = reshape(im,[hd(2) hd(3) hd(4)]); 
    else
    im = reshape(im,[hd(2) hd(3) 2 hd(4)]); 
    im = squeeze(complex(im(:,:,1,:),im(:,:,2,:)));
end



end