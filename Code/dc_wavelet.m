function dcData = dc_wavelet(dcfile)

    [m,n] = size(dcfile); % 4096 x 80
    pxl = sqrt(m);
    nw = m/4; % wavelet resolution
    dcData = zeros(nw,n);
    
    for k = 1:n
        X = im2double(reshape(dcfile(:,k),pxl,pxl));
        [~,cH,cV,~]=dwt2(X,'haar');
        cod_cH1 = rescale(abs(cH));
        cod_cV1 = rescale(abs(cV));
        cod_edge = cod_cH1+cod_cV1;
        dcData(:,k) = reshape(cod_edge,nw,1);
    end
end
