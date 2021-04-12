function Lrgb_final=pre_process1(provoli_eikonwn,I,rad,sirikn_matrix,pix_num)


hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate'); %replicate για να μην δημιουργεί η συνέλιξη με το φίλτρο μαύρες γραμμές
Ix = imfilter(double(I), hx, 'replicate'); %στα άκρα, γιατί θεωρεί ότι εκτός εικόνας τα πίξελ έχουν τιμή 0, το replicate
gradmag = sqrt(Ix.^2 + Iy.^2);             %αντιγράφει τα άκρα της εικόνας όπως είναι. Υπολογίζω την παράγωγο και βρίσκω το πλάτος της

% Για να βρω τους markers, χρησιμοποιώ opening-by-reconstruction και
% closing-by-reconstruction, που θα δημιουργήσει επίπεδα μέγιστα μέσα σε
% κάθε αντικείμενο ( εντοπίζεται με imregionalmax).
% Opening είναι ένα erosion που ακολουθείται από dilation,
% ενώ opening-by-reconstruction είναι ένα erosion που ακολουθείται morphological reconstruction. 


%rad=5; %ακτίνα δίσκου αντικειμένων που θέλω να αφαιρέσω με το opening
se = strel('disk', rad); 


Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);


Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);


%Υπολογισμός του regional maxima για να φτιάξω foreground markers.

Iobrcbr1=Iobrcbr;
Iobrcbr1(find(Iobrcbr1<30/255))=0;
fgm = imregionalmax(Iobrcbr1);


I2 = I;
I2(fgm) = 255;


% Bλέπω οτι κάποιοι foreground markers εκτείνονται μέχρι τις άκρες των αντικειμένων
% Για αυτό πρέπει να καθαρίσω τις ακμές από τις περιοχές των markers και να
% τους συρρικνωσω. Εφαρμόζω closing και έπειτα erosion.

se2 = strel(sirikn_matrix);
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);

% Eξαλείφω περιοχές με λιγότερα απο 5 πιξελ
fgm4 = bwareaopen(fgm3,pix_num);
I3 = I;
I3(fgm4) = 255;

if provoli_eikonwn
figure;
subplot(1,2,1);
imshow(I);
title('Original image (I)');
subplot(1,2,2);
imshow(I2);
title('Regional maxima superimposed on original image (I2)');
end;


thresh=graythresh(Iobrcbr);
%thresh=triangle_th(imhist(Iobrcbr),256);


bw = im2bw(Iobrcbr,thresh);
%σκέψου να χρησιμοποιήσεις τη συνάρτηση τριγώνου για thresholding στους
%μεγάλους πυρήνες


if provoli_eikonwn
figure;
subplot(2,2,1);
imshow(Iobr);
title('opening-by-reconstruction (Iobr)');
subplot(2,2,2);
imshow(Iobrcbr);
title('Opening-closing by reconstruction (Iobrcbr)');
subplot(2,2,3);
imshow(Iobrcbr1);
title('Iobrcbr with pixel under 30 insensity being total black');
subplot(2,2,4);
imshow(bw);
title('Thresholded Iobrcbr (bw)');
end;

%Τα background pixels είναι μάυρα, όμως δεν θέλω οι background markers να είναι πολύ κοντά στις ακμές των 
%αντικειμένων που θέλω να χωρίσω. Μπορώ να λεπτύνω το background υπολογίζοντας το σκελετό του foreground.

D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;

%modify the gradient magnitude image so that its only regional minima occur at foreground and background marker pixels.
gradmag2 = imimposemin(gradmag, bgm | fgm4);

L = (watershed(gradmag2));

if provoli_eikonwn
figure;
subplot(2,2,1);
imshow(gradmag2);
title('Modifing gradient');
subplot(2,2,2);
imshow(bgm);
title('Watershed ridge lines (bgm)');
subplot(2,2,3);
imshow(fgm4);
title('deleting areas with less pixels than a given value'); 
subplot(2,2,4);
imshow(L);
title('watershed segmentation');
end;

I4 = I;
I4(imdilate(L == 0, sirikn_matrix) | bgm | fgm4) = 255;



Lrgb = label2rgb(L, 'white', 'k', 'shuffle');


Lrgb2=rgb2gray(Lrgb);
Lrgb2=~Lrgb2;
Lrgb2_filled=imfill(Lrgb2,'holes');


Lrgb_final=Lrgb2_filled-Lrgb2;

if provoli_eikonwn
figure;
subplot(2,2,1);
imshow(I4);
title('Markers and object boundaries superimposed on original image (I4)');
subplot(2,2,2);
imshow(Lrgb2);
title('boundaries of the image');
subplot(2,2,3);
imshow(Lrgb2_filled);
title('image with filled holes');
subplot(2,2,4);
imshow(Lrgb_final);
title('final image as filled image-boundaries');
end;



end