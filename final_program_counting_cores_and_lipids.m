%% GUI
%επιλογή φακέλου αποθήκευσης αποτελεσμάτων
saving_path=uigetdir('C:\Users\','Choose folder to save the results');

%επιλογή τρόπου εκτέλεσης
choice = questdlg('Choose Procedur','Execution Procedure','One image','Folder',...
    'Folder with subfolders','Folder with subfolders');
% Handle response
switch choice
    case 'One image'
       ektelesh_mia_eikona = 1;
       ektelesh_fakelo=0;
       ektelesh_pollwn_fakelwn=0;
    case 'Folder'
       ektelesh_mia_eikona=0;
       ektelesh_fakelo=1;
       ektelesh_pollwn_fakelwn=0;
    case 'Folder with subfolders'
       ektelesh_mia_eikona=0;
       ektelesh_fakelo=0;
       ektelesh_pollwn_fakelwn=1;
end


%για επιλογή εκτέλεσης μίας εικόνας, επιλογή της συγκεκριμένης εικόνας
if ektelesh_mia_eikona
    [filename_lipids, pathname_lipids]=uigetfile({'*.jpg;*.tif;*.png;*','All Image Files'},...
    'Choose image with lipids','C:');

    [filename_cores, pathname_cores]=uigetfile({'*.jpg;*.tif;*.png;*','All Image Files'},...
    'Choose image with cores','C:');  
    
%επιλογή εμφάνισης ή μη των αποτελεσμάτων
    choice = questdlg('Display images of processing steps?','Images of processing steps','All images','Only final images','None','None');
    switch choice
       case 'All images'
        provoli_eikonwn = 1;
        provoli_telikwn_eikonwn=1;
       case 'Only final images'
        provoli_eikonwn=0;
        provoli_telikwn_eikonwn=1;
       case 'None'
        provoli_eikonwn = 0;
        provoli_telikwn_eikonwn=0;
    end
else
    %επιλογή φορμάτ εικόνων για εκτέλεση πολλών εικόνων
    choice = questdlg('Choose format of images;','Choose format','jpg','bmp','tiff','tiff');
    switch choice
      case 'jpg'
       format_jpg=1;
       format_bmp=0;
       format_tiff=0;
      case 'bmp'
       format_jpg=0;
       format_bmp=1;
       format_tiff=0;
      case 'tiff'
       format_jpg=0;
       format_bmp=0;
       format_tiff=1;
    end
    if ektelesh_fakelo
      path_name=uigetdir;
      if format_jpg   
       images_list = dir(fullfile(path_name,'*.jpg'));
      elseif format_bmp
       images_list = dir(fullfile(path_name,'*.bmp'));
      elseif format_tiff
       images_list = dir(fullfile(path_name,'*.TIFF'));
      end
      names = {images_list.name}';

      c1=0;
      c2=0;
      for i=1:length(names)
        if isempty(strfind(names{i},'NR_'))
         c1=c1+1;
         filename_cores_cell{c1,1}=names{i};
         split{c1}=strsplit(filename_cores_cell{c1},{'H_','.'});
         split1(c1)=split{c1}(2);
         split2(c1)=split{c1}(1);
         split_array(c1)=str2num(split1{c1});
         filename_cores_cell{c1,2}=split_array(c1);
       else
         c2=c2+1;
         filename_lipids_cell{c2,1}=names{i}; 
         splitl{c2}=strsplit(filename_lipids_cell{c2},{'NR_','.'});
         split1l(c2)=splitl{c2}(2);
         split_arrayl(c2)=str2num(split1l{c2});
         filename_lipids_cell{c2,2}=split_arrayl(c2);
       end;
      end;
 
      filename_cores_cell=sortrows(filename_cores_cell,2);
      filename_lipids_cell=sortrows(filename_lipids_cell,2);
      for i=1:length(split_array)
        filename_cores{i}=filename_cores_cell{i,1};
        filename_lipids{i}=filename_lipids_cell{i,1};
      end;    
    end;
    if ektelesh_pollwn_fakelwn
      path_name=uigetdir;
      all_folders=genpath(path_name);
      path_sep=strfind(all_folders,';');
      c1=0;
      c2=0;
      for i=1:length(path_sep)-1
        path_sub{i}=all_folders(path_sep(i)+1:path_sep(i+1)-1);
        if  strfind(path_sub{i},'DAPI')
          c1=c1+1;
          filename_cores{c1}=path_sub{i};
        end;
        if strfind(path_sub{i},'RFP')
          c2=c2+1;
          filename_lipids{c2}=path_sub{i};
        end;
      end;    
    end;
    choice = questdlg('Display images of processing steps?','Images of processing steps','All images in every loop',...
      'Only final images in the end of the program','None','None');
        % Handle response
    switch choice
       case 'All images in every loop'
       provoli_eikonwn = 1;
       provoli_telikwn_eikonwn=1;
       case 'Only final images in the end of the program'
       provoli_eikonwn=0;
       provoli_telikwn_eikonwn=1;
       case 'None'
       provoli_eikonwn = 0;
       provoli_telikwn_eikonwn=0;
    end
end;
 choice = questdlg('Kind of Cores''s images','Cores''s images','Small Cores',...
     'Big Cores','Big Cores');
% Handle response
switch choice
    case 'Small Cores'
       small_cores = 1;
    case 'Big Cores'
      small_cores=0;
end
 
%εισαγωγή παραμέτρων από το χρήση ή επιλογή default τιμών
w= [1 1 1; 1 -8 1; 1 1 1]; 
   rad1=3; %ακτίνα δίσκου αντικειμένων που θέλω να αφαιρέσω με το opening 
   sirikn_matrix1=ones(1,1);
   pix_num1=5;
   rad2=5; %ακτίνα δίσκου αντικειμένων που θέλω να αφαιρέσω με το opening 
   sirikn_matrix2=ones(3,3);
   pix_num2=5;
     rad4=15; %ακτίνα δίσκου αντικειμένων που θέλω να αφαιρέσω με το opening 
  min_bright=200;
  t1=0.225;
  t2=0.825;
  t3=0.85;
  t4=0.875;
  t5=0.9;
  t6=0.925;
  n1=1;
  n2=5;
  n3=10;
  n4=15;
  n5=20;
  ec_level=0.7;
  round_level=0.7;
  circ_level=0.7;
choice = questdlg('Choose default parameters;','Parameters','Yes','No','No');
% Handle response
switch choice
    case 'Yes'
      default_times=1;
    case 'No'
      default_times=0;
end

if default_times
   max_area=250;
   min_lipid=5;
   max_lipid=1500;
else
    
    
    prompt = {'max area of small lipids during processing (default=250,suggested=150-350):',...
        'min area of lipids after processing (default=5,suggested=4-10):',...
        'max area of lipids after processing (default=1500,suggested=1000-2000):'};
dlg_title = 'Processing Image with Lipids'; 
answer = inputdlg(prompt,dlg_title,[1 50;1 50;1 50]);

  if isempty(answer{1})
  max_area=250; 
  else
  max_area=str2num(answer{1});
  end; 
  if isempty(answer{2})
  min_lipid=5; 
  else
  min_lipid=str2num(answer{2});
  end;
  if isempty(answer{3})
  max_lipid=1500;
  else
  max_lipid=str2num(answer{3});
  end;

end;

%%



error_counter=0;
if class(filename_lipids)=='cell'
loops=numel(filename_lipids);
else
loops=1;
end;
for k=1:loops
    tic
    disp(sprintf('<< Processing Image %d >>',k));
    try
       
       
       if ektelesh_mia_eikona
        image_red1=imread(strcat(pathname_lipids,filename_lipids));
        image_blue1=imread(strcat(pathname_cores,filename_cores));
        if size(image_red1,3)==4
            image_red1=image_red1(:,:,1:3);
            image_blue1=image_blue1(:,:,1:3);
        end;
       elseif ektelesh_fakelo
        image_blue1 = imread( fullfile( path_name, filename_cores{k} ) );
        image_red1 = imread( fullfile( path_name, filename_lipids{k} ) );
          if size(image_red1,3)==4
            image_red1=image_red1(:,:,1:3);
            image_blue1=image_blue1(:,:,1:3);
          end;
       elseif ektelesh_pollwn_fakelwn
         if format_tiff
          image_red1=imread(fullfile(filename_lipids{k},'00000.TIFF')); 
          image_blue1=imread(fullfile(filename_cores{k},'00000.TIFF')); 
          image_red1=image_red1(:,:,1:3);
          image_blue1=image_blue1(:,:,1:3);
         end
         if format_bmp
          image_red1=imread(fullfile(filename_lipids{k},'00000.bmp')); 
          image_blue1=imread(fullfile(filename_cores{k},'00000.bmp')); 
         end
         if format_jpg
          image_red1=imread(fullfile(filename_lipids{k},'00000.jpg')); 
          image_blue1=imread(fullfile(filename_cores{k},'00000.jpg')); 
         end
       end;
       
       image_red=image_red1(:,:,1);
       image_blue=image_blue1(:,:,3);

  

%% LAPLACIAN FILTER

%w= [1 1 1; 1 -8 1; 1 1 1];

lap1_red_w=imfilter(image_red,w,'replicate');
image_red2=im2double(image_red);
lap2_red_w=imfilter(image_red2,w,'replicate');
lap_red_final_w=image_red2-lap2_red_w;



lap1_blue_w=imfilter(image_blue,w,'replicate');
image_blue2=im2double(image_blue);
lap2_blue_w=imfilter(image_blue2,w,'replicate');
lap_blue_final_w=image_blue2-lap2_blue_w;


%% MEDIAN FILTER FOR NOISE REDUCTION

image_red_med=medfilt2(lap_red_final_w,'symmetric');
image_blue_med=medfilt2(lap_blue_final_w,'symmetric');

if provoli_eikonwn
     
    image_red_plot=image_red1(:,:,1:3);
    image_blue_plot=image_blue1(:,:,1:3);
    
    lap_red_final_w_plot(:,:,1)=lap_red_final_w;
    lap_red_final_w_plot(:,:,2)=0*lap_red_final_w;
    lap_red_final_w_plot(:,:,3)=0*lap_red_final_w;
    
    lap_blue_final_w_plot(:,:,1)=0*lap_blue_final_w;
    lap_blue_final_w_plot(:,:,2)=0*lap_blue_final_w;
    lap_blue_final_w_plot(:,:,3)=lap_blue_final_w;
    
    image_red_med_plot(:,:,1)=image_red_med;
    image_red_med_plot(:,:,2)=0*image_red_med;
    image_red_med_plot(:,:,3)=0*image_red_med;
    
    image_blue_med_plot(:,:,1)=0*image_blue_med;
    image_blue_med_plot(:,:,2)=0*image_blue_med;
    image_blue_med_plot(:,:,3)=image_blue_med;

figure;
subplot(2,3,1);
imshow(image_red_plot);
title('red image');
subplot(2,3,2);
imshow(lap_red_final_w_plot);
title('red image enhanced by w laplacian filtering');
subplot(2,3,3);
imshow(image_red_med_plot);
title('symmetric median filtering ');
subplot(2,3,4);
imshow(image_blue_plot);
title('blue image');
subplot(2,3,5);
imshow(lap_blue_final_w_plot);
title('blue image enhanced by w laplacian filtering');
subplot(2,3,6);
imshow(image_blue_med_plot);
title('symmetric median filtering ');
end;





%% LIPIDS

I1=image_red_med; %max_area=300;
Lrgb_final_red_small_lipids=pre_process1(provoli_eikonwn,I1,rad1,sirikn_matrix1,pix_num1);
Lrgb_final_red_small_lipids_labeled=bwlabel(Lrgb_final_red_small_lipids);

stats1=regionprops(Lrgb_final_red_small_lipids_labeled,'Area');

for i=1:length(stats1);
    if stats1(i).Area>max_area
       Lrgb_final_red_small_lipids_labeled(find(Lrgb_final_red_small_lipids_labeled==i))=0;
    end;
end;

%περιοχές με αριθμό πίξελ τόσο, θα εξαλειφθούν
I2=image_red-image_blue;%αφαίρεση από την κόκκινη την μπλε εικόνα που έχει κοινό θόρυβο
Lrgb_final_red_big_lipids=pre_process1(provoli_eikonwn,I2,rad2,sirikn_matrix2,pix_num2); 
Lrgb_final_red_big_lipids=imerode(Lrgb_final_red_big_lipids,strel('disk',4));

Lrgb_total=Lrgb_final_red_big_lipids|Lrgb_final_red_small_lipids_labeled;

element=strel('disk',2);
Lrgb_total=imerode(Lrgb_total,strel('disk',1));
Lrgb_total=bwlabel(Lrgb_total);
stats_lipids_cell{k}=regionprops(Lrgb_total,'Area');
 
counter=0;
for i=1:length(stats_lipids_cell{k});
    if cell2mat(struct2cell(stats_lipids_cell{k}(i)))>max_lipid||cell2mat(struct2cell(stats_lipids_cell{k}(i)))<=min_lipid
       Lrgb_total(find(Lrgb_total==i))=0;
       counter=counter+1;
       positions(counter)=i;
    end;
end;

if counter>0
stats_lipids_cell{k}(positions)=[];
end;

    num_of_lipids(k)=max(max(bwlabel(Lrgb_total)));
    stats_lipids_array=sort(cell2mat(struct2cell(stats_lipids_cell{k})));
    mean_area_lipids(k)=mean(stats_lipids_array);
    median_area(k)=median(stats_lipids_array);
    sd(k)=std(stats_lipids_array); 
    total_area_lipids(k)=sum(stats_lipids_array);
    clear stats_lipids_array positions



I_rgb_lipids(:,:,1) = 255*Lrgb_total;
I_rgb_lipids(:,:,2) = 0*Lrgb_total;
I_rgb_lipids(:,:,3) = 0*Lrgb_total;

I_rgb_lipids_green(:,:,1) = 0*Lrgb_total;
I_rgb_lipids_green(:,:,2) = 255*Lrgb_total;
I_rgb_lipids_green(:,:,3) = 0*Lrgb_total;

%% CORES 

if small_cores


rad5=10;
sirikn_matrix5=ones(5,5);
pix_num5=150;

se=strel('disk',rad4);
Ie = imerode(image_red_med, se);
Iobr = imreconstruct(Ie,image_red_med);
Ic=pre_process1(0,image_blue_med-Iobr,rad5,sirikn_matrix5,pix_num5);

Ic=bwlabel(Ic);

stats5=regionprops(Ic);
Ic1=Ic;
for i=1:length(stats5);
    if stats5(i).Area>10000||stats5(i).Area<100
       Ic1(find(Ic1==i))=0;
    end;
end;

if provoli_eikonwn
figure;
subplot(1,2,1);
imshow(Ic);
subplot(1,2,2);
imshow(Ic1);
end

num_of_cores_var_th(k)=max(max(bwlabel(Ic1)));
 
I_rgb_cores(:,:,1) = 0*Ic1;
I_rgb_cores(:,:,2) = 0*Ic1;
I_rgb_cores(:,:,3) = Ic1;

I_rgb_cores_cyan(:,:,1) = 0*Ic1;
I_rgb_cores_cyan(:,:,2) = Ic1;
I_rgb_cores_cyan(:,:,3) = Ic1;
I_rgb_cores_cyan=im2double(I_rgb_cores_cyan);

    
    
else
se=strel('disk',rad4);
Ie = imerode(image_red_med, se);
Iobr = imreconstruct(Ie,image_red_med);


image_blue_minus_red=image_blue_med-Iobr;
    

count_brightness(k)=100*sum(sum(image_blue_minus_red>=min_bright/255))/(size(image_blue_minus_red,1)*size(image_blue_minus_red,2));

n01=0.0001;
n02=0.0005;
n03=0.001;
n04=0.007;
n05=0.01;
if count_brightness(k)<=n01
level=0.0005;   
elseif count_brightness(k)>n01&&count_brightness(k)<=n02
level=0.0001;
elseif count_brightness(k)>n02&&count_brightness(k)<=n03
level=0.005;
elseif count_brightness(k)>n03&&count_brightness(k)<=n04
level=0.1; 
elseif count_brightness(k)>n04&&count_brightness(k)<=n05
level=0.14;
elseif count_brightness(k)>n05&&count_brightness(k)<=n1    
level=t1;
elseif count_brightness(k)>n1&&count_brightness(k)<=n2
level=t2;
elseif count_brightness(k)>n2&&count_brightness(k)<=n3
level=t3;
elseif count_brightness(k)>n3&&count_brightness(k)<=n4
level=t4;
elseif count_brightness(k)>n4&&count_brightness(k)<=n5
level=t5;
else
level=t6;
end;


image_blue_thresholded=im2bw(image_blue_minus_red,level);

Ic2=bwlabel(image_blue_thresholded);
stats6=regionprops(Ic2);
for i=1:length(stats6);
    if stats6(i).Area<20
       Ic2(find(Ic2==i))=0;
    end;
end;

clear stats6

image_blue_dilated=imdilate(Ic2,strel('disk',4));
image_blue_filled=imfill(image_blue_dilated,'holes');


if provoli_eikonwn
figure;
subplot(2,3,1);
imshow(image_blue_med);
title('image blue');
subplot(2,3,2);
imshow(Iobr);
title('image red after reconstruction opening with disk 4');
subplot(2,3,3);
imshow(image_blue_minus_red);
title('blue image-red reconstructed image');
subplot(2,3,4);
imshow(image_blue_thresholded);
title('thresholded image blue');
subplot(2,3,5);
imshow(image_blue_dilated);
title('dilated thresholded image blue');
subplot(2,3,6);
imshow(image_blue_filled);
title('filled blue image');
end;

Ie1 = imerode(image_blue_filled, se);
Iobr2 = imreconstruct(Ie1,image_blue_filled);


stats2=regionprops(Iobr2,'Area','Centroid');


%απαλοιφή μεγάλων περιοχών
Iobr1=Iobr2;
counter_bc=0;
for i=1:length(stats2);
    if stats2(i).Area>10000
       counter_bc=counter_bc+1;
       Iobr1(find(Iobr2==i))=0;
       positions_bc(counter_bc)=i;
    end;
end;

if counter_bc>0
stats2(positions_bc)=[];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:length(stats2)
  center1(i,:)=stats2(i).Centroid;
  area(i,1)=stats2(i).Area;
  rad3(i)=sqrt(area(i)/pi);
end;
  
mask = createCirclesMask(Iobr1,center1,rad3);

I=Iobr1&mask;


I_rgb_cores(:,:,1) = 0*I;
I_rgb_cores(:,:,2) = 0*I;
I_rgb_cores(:,:,3) = I;

I_rgb_cores_cyan(:,:,1) = 0*I;
I_rgb_cores_cyan(:,:,2) = I;
I_rgb_cores_cyan(:,:,3) = I;
I_rgb_cores_cyan=im2double(I_rgb_cores_cyan);


stats3=regionprops(I,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','EquivDiameter');

  
  boundary=bwboundaries(I,'noholes');
  for i=1:length(boundary)
  center2(i,:)=stats3(i).Centroid;
  area(i,1)=stats3(i).Area;
  delta_sq = diff(boundary{i}).^2;
  perimeter(i,1) = sum(sqrt(sum(delta_sq,2)));
  roundness(i,1)=4*(area(i))/(pi*stats3(i).MajorAxisLength^2);
  circularity(i,1)= 4*pi*area(i)/perimeter(i)^2;
  eccentricity(i,1)=1-stats3(i).Eccentricity;% ορίζω εκκεντρικότητα 1-εκκεντρικότητα
  end;
  


  %Roundness=(4*Area)/(π*(max_diameter)^2)
  %Circularity=(4*π*Αrea)/perimeter^2
  %Eccentricity= ratio of the distance between the foci of the ellipse,that has 
  %the same second-moments as the region, and its major axis length.
  %Area with eccentricity 0 is a circle, while area with eccentricity 1 is
  %a line
  %Εγώ την ορίζω 1-εκκεντρικότητα ώστε να είναι το 1 κύκλος και το 0 γραμμή
  %για να είναι συγκρίσιμα με το roundness και το circularity
  
ec_level=0.7;
round_level=0.7;
circ_level=0.7;
num_of_cores_var_th(k)=0;
for i=1:length(stats3);
if (eccentricity(i)<ec_level)&&(roundness(i)<round_level)&&(circularity(i)<circ_level)
num_of_cores_var_th(k)=num_of_cores_var_th(k)+2;    
else
num_of_cores_var_th(k)=num_of_cores_var_th(k)+1;
end;
end;

  
end%απο το If επιλογης μικρών πυρήνων

percentage(k)=100*num_of_lipids(k)/num_of_cores_var_th(k);


if provoli_telikwn_eikonwn&&provoli_eikonwn
figure;
subplot(1,2,1);
imshow(image_red1);
title('Lipids');
subplot(1,2,2);
imshow(image_red1);
hold on
himage = imshow(I_rgb_lipids_green);
himage.AlphaData = 0.5;
title('Lipids after processing(green)');
hold off;

figure;
subplot(1,2,1);
imshow(image_blue1);
title('Cores');
subplot(1,2,2);
imshow(image_blue1);
hold on
himage = imshow(I_rgb_cores_cyan);
himage.AlphaData = 0.5;
title('Cores after processing(cyan)');
hold off;
end;

if provoli_telikwn_eikonwn
figure;
subplot(1,2,1);
imshow(image_blue1+image_red1);
subplot(1,2,2);
imshow(I_rgb_cores+I_rgb_lipids);
visib_fig='on';
else
visib_fig='off';    
end;

[counts,binCenters]=hist(cell2mat(struct2cell(stats_lipids_cell{k})), 2000);
h=figure;
set(h,'Visible',visib_fig);
subplot(2,1,1);
bar(binCenters,counts);
hold on
bar(mean_area_lipids(k),0.8*max(counts),'r','linewidth',2,'EdgeColor','r');
bar(median_area(k),0.8*max(counts),'g','linewidth',2,'EdgeColor','g');
bar(sd(k),0.8*max(counts),'m','linewidth',2,'EdgeColor','m');
legend({'Number of Lipids','Mean Value','Median Value','Standard Deviation'},'Location','northeast');
title('Lipids Area Distribution')
xlabel('Area (pixels)');
ylabel('Frequency');
subplot(2,1,2);
bar(log10(binCenters),counts);
title('log(Lipids Area) Distribution');
xlabel('log(Area)');
ylabel('Frequency');


%% SAVING RESULTS
imwrite(image_blue1+image_red1,fullfile(saving_path,sprintf('%d_original_image.jpg',k)));
imwrite(I_rgb_cores+I_rgb_lipids,fullfile(saving_path,sprintf('%d_final_image.jpg',k)));
saveas(h,fullfile(saving_path,sprintf('histogram_image_%d.jpg',k)));
%%


error_position(k)=0;
timer_loop(k)=toc;

if timer_loop(k)<60
disp(sprintf('<< Elapsed Time %d seconds >>',round(timer_loop(k))));
else
disp(sprintf('<< Elapsed Time %s minutes >>',char(vpa(timer_loop(k)/60,3))));
end;


%toc;
if ektelesh_mia_eikona
    break;
else
  if  provoli_eikonwn
  promptMessage = sprintf('Do you want to Continue processing,\nor Cancel to abort processing?');
  button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
  if strcmpi(button, 'Cancel');
  return; % Or break or continue
  end
  end;
end;

clear counts binCenters stats1 stats2 positions_bc stats3 stats5 eccentricity circularity area rad3 boundary center1 center2 perimeter roundness delta_sq
  

    catch ME
    timer_loop(k)=toc;    
    disp(sprintf('<< Error in Image %d >>',k));
    disp(sprintf('<< Elapsed Time %d seconds >>',round(timer_loop(k))));
    error_position(k)=k;
    end;

end; %end του k loop


k_matrix=1:k;
if error_position==k_matrix
    if k>1
    disp('<< Error in all Images >>');
    end;
return;
end;

for j=find(error_position>0)
num_of_lipids(j)=0;
num_of_cores_var_th(j)=0;
total_area_lipids(j)=0;
percentage(j)=0;
mean_area_lipids(j)=0;
median_area(j)=0;
sd(j)=0;
end;

if ~ektelesh_mia_eikona
h2=figure;
set(h2,'Visible','off');
subplot(3,1,1);
bar((1:k),[num_of_lipids;num_of_cores_var_th; total_area_lipids/100]');%,0.6)%,'r');
ylabel('Number');
xlabel('Number of Image');
legend({'Lipids','Cores','Area of Lipids/100'},'Location','north','Orientation','vertical','FontSize',7);

subplot(3,1,2);
bar((1:k),[percentage;(total_area_lipids)./num_of_cores_var_th ; 100000.*(num_of_lipids./total_area_lipids) ]');%,0.25);%,'g');
ylabel('100%');
xlabel('Number of Image');
legend({'Lipids/Cores','Lipids(Area)/Cores','1000*Lipids/Area of Lipids'},'Location','north','Orientation','vertical','FontSize',7);

subplot(3,1,3);
bar((1:k),[mean_area_lipids;median_area;sd]');
ylabel('Area of Lipids (pixels)');
xlabel('Number of Image');
hold on;
legend({'Mean Value','Median Value','Standard Deviation'},'Location','north','Orientation','vertical','FontSize',7);
saveas(h2,fullfile(saving_path,'barplots2.jpg'));

end;


for i=find(k_matrix~=error_position)  %length(num_of_lipids)
if ektelesh_fakelo 
file{i,1}=cell2mat(split2(i));
else
file{i,1}=sprintf('%s_%d','image',i);
end;
file{i,2}=num_of_lipids(i);
file{i,3}=num_of_cores_var_th(i);
file{i,4}=percentage(i);
file{i,5}=total_area_lipids(i)/num_of_cores_var_th(i);
file{i,6}=total_area_lipids(i);
file{i,7}=100*num_of_lipids(i)/total_area_lipids(i);
file{i,8}=mean_area_lipids(i);
file{i,9}=median_area(i);
file{i,10}=sd(i);
end;

  if i==k
    i=k;
  else
    i=error_position(end);
  end;


num_of_lipids(error_position>0)=[];
num_of_cores_var_th(error_position>0)=[];
percentage(error_position>0)=[];
total_area_lipids(error_position>0)=[];
mean_area_lipids(error_position>0)=[];
median_area(error_position>0)=[];
sd(error_position>0)=[];

file{i+2,1}='Mean_Value';
file{i+3,1}='Median_Value';
file{i+4,1}='Standard_Deviation';
file{i+2,2}=mean(num_of_lipids);
file{i+3,2}=median(num_of_lipids);
file{i+4,2}=std(num_of_lipids);
file{i+2,3}=mean(num_of_cores_var_th);
file{i+3,3}=median(num_of_cores_var_th);
file{i+4,3}=std(num_of_cores_var_th);
file{i+2,4}=mean(percentage);
file{i+3,4}=median(percentage);
file{i+4,4}=std(percentage);
file{i+2,5}=mean(total_area_lipids./num_of_cores_var_th);
file{i+3,5}=median(total_area_lipids./num_of_cores_var_th);
file{i+4,5}=std(total_area_lipids./num_of_cores_var_th);
file{i+2,6}=mean(total_area_lipids);
file{i+3,6}=median(total_area_lipids);
file{i+4,6}=std(total_area_lipids);
file{i+2,7}=mean(100*num_of_lipids./total_area_lipids);
file{i+3,7}=median(100.*num_of_lipids./total_area_lipids);
file{i+4,7}=std(100.*num_of_lipids./total_area_lipids);

for i=find(k_matrix==error_position) 
if ektelesh_fakelo 
file{i,1}=cell2mat(split2(i));
else
file{i,1}=sprintf('%s_%d','image',i);
end;  
file{i,2}='ERROR';
file{i,3}='ERROR';
file{i,4}='ERROR';
file{i,5}='ERROR';
file{i,6}='ERROR';
file{i,7}='ERROR';
file{i,8}='ERROR';
file{i,9}='ERROR';
file{i,10}='ERROR';
end;


T=cell2table(file,'VariableNames',{'Image' 'Lipids' 'Cores' 'Percentage_Lipids_Cores_x100'...
    'Total_Area_Num_of_Cores' 'Total_Area_Lipids' 'Lipids_Num_per_Area' 'Mean_Area_Lipids'...
    'Median_Area_Lipids' 'StD_Area_Lipids' });
writetable(T,fullfile(saving_path,'results.xls'));

if sum(timer_loop)<60
disp(sprintf('<< Programm Completed in %d seconds >>',round(sum(timer_loop))));
else 
disp(sprintf('<< Programm Completed in %s minutes >>',char(vpa(sum(timer_loop)/60,3))));
end;
