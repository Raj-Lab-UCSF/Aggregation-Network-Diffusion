function display_img(q, disp_text, thr_ptile)
% by A Raj: displays 2d image as single figure, 3D as a collage
% Usage: display_img(img, disp_text)

SIZ = size(q);
if nargin<2 || isempty(disp_text)
    disp_text = 'image';
end
if nargin<3 || isempty(thr_ptile)
    thr_ptile = [0.01, 0.99];
end

% scale images to saturate at 99%tile
sq = abs(q(:));
sq = sort(sq, 'Ascend');
thr_lo = sq(max(1,round(length(sq)*thr_ptile(1))));
sq = sq(sq>=thr_lo);
% sq = abs(q(:));
% sq = sq(sq>0);
sq = sort(sq, 'Ascend');
thr_lo = sq(max(1,round(length(sq)*thr_ptile(1))));
thr_hi = sq(min(length(sq), round(length(sq)*thr_ptile(2))));
if thr_hi > thr_lo
    q(q>thr_hi) = thr_hi;
    q(q<thr_lo) = thr_lo;
end
if length(SIZ)==2
    figure('Name', disp_text);  imagesc(q);  %colormap(gray);
    blowup = max(SIZ(1), 512)/SIZ(1);
    title(disp_text);
    truesize(round(blowup*SIZ));
elseif length(SIZ)==3
    [nrows, ncols, nslices] = size(q);
    blowup = max([SIZ(1), SIZ(2), 256])/max(SIZ(1), SIZ(2));
    if nslices<7
        dd = 1:nslices;
    elseif nslices==7
        dd = 1:6;
    elseif nslices>7
        dd = round(linspace(2, nslices-1, 6));
        dd = unique(dd);
    end
    dr = sqrt(length(dd));
    dr = floor(dr);
    dc = ceil(length(dd)/dr);
    figure('Name', disp_text);  %colormap(gray);
    for i = 1:length(dd)
         subplot(dr, dc, i); imagesc(q(:,:,dd(i)));  %truesize(round(blowup*SIZ(1:2))); %title(['slice # ' num2str(i)]);
         %subplot(dr, dc, i); imshow(q(:,:,dd(i)),'InitialMagnification','fit');  
         title( ['slice #' num2str(dd(i))] ); 
         if thr_hi > thr_lo, caxis([thr_lo thr_hi]); end
         colormap jet;
         axis off;
    end
elseif length(SIZ)==4
   disp('Read 4D block - will not display');
end
