%This is the code for Martha Bhattacharya's MCB585 project
%The purpose of this script is to create regions around neurons in an
%imaging data set, and extract their intensity over the image set (time course). 
%It makes use of segmentation and loops.

%In most cases I have commented out the creation of figures to mininimze
%pop-ups. The last figure with adjusted regions will still appear.

%First I will open up the last image in an image series, which is the image that I
%want to use for creating regions.

%Need to do this by capturing names, number of tiffs, and
%finding the last tiff in numerical order.

%Here I will try to create a dialog box to help put in the appropriate
%folder name for reading tiff files
prompt = {'What folder are the tif images in? Use full path:'};
dlgtitle = 'File Location';
dims = [1 50];
definput = {'O:\MBdata\Ca_imaging'};
FileLoc = inputdlg(prompt,dlgtitle,dims,definput);

%convert to a character output
FileLoc = char(FileLoc);

%make list of file names; create an empty struct for data for later
ImageDir = dir(FileLoc);
stats = struct;
NumFiles = size(ImageDir,1);
ImageDir = ImageDir(3:NumFiles,:);
%note that this gives output of a few extra files (.,..,a folder,one other)

%Now need to make the first "worked on" image the last tiff in the folder.
%I am going to try to write a function that counts tiffs
numtiffs = 0;
for i = 1:size(ImageDir,1)
    ImName = ImageDir(i).name;
    istif = regexp(ImName,'.tif');
    if (istif)>0
        numtiffs = numtiffs + 1;
    end
end
%j = 1:size(ImageDir,1)
for j = 1:numtiffs
    ImName = ImageDir(j).name;
    istif = regexp(ImName,'.tif');
    if (istif)>0
        frame = regexp(ImName,'.+_t(\d{2}).+.tif','tokens','once');
        %frame = regexp(ImName,'.+_t(\d{3}).+','tokens','once');
    end
    if str2num(frame{1}) == numtiffs
       HiKframe = ImName;
    end
end

%This gives a character vector that is the name of the file.
%Now open the file with this name
fullname = [FileLoc '\' HiKframe];
im = imread(fullname);

%Now that we have the correct image open, we will do the segmentation on
%the green image. Ideally we'd do it on the ratio, but when I divide green
%by red, it does the division but rounds to the nearest integer, losing
%data. I can't find the command to NOT round. So will just do green for
%now.

     %split into red, green, blue images (Fura-2 imaging only!)
        %[R,G,B] = imsplit(im);
     %the below command also was changed from G to im for Fluo-4 single
     %color imaging
        imT = RollingBall(im,5,50);
        %imT = RollingBall(im(:,:,3),5,50);
        imTG = imgaussfilt(imT,1);
        Thresh = im2double(multithresh(imTG)./10); 
        %Otsu's method threshold is too high, lose most cells, so divide by 10 instead
        %will have to test whether this is an okay threshold for majority of images
        BI = imbinarize(imTG,Thresh);
        nopixels = medfilt2(BI,[12 12]);
        %I had to optimize the dimensions in the median filter so that it got rid
        %of just the smallest non-cellular pixels
        %figure, imshowpair(BI,nopixels,'montage');
        %This looks fairly good.
        
        %BI = imclearborder(BI);
        
        BI = SplitClumps(BI,1);
        [L,num] = bwlabel(BI); %L makes each unique clump a number or label
        figure, imshow(label2rgb(L, 'jet','k','shuffle'));
       
%Great! Now there are still some large clumps that were not broken, and
%some very small regions that are not neurons. Still need to exclude, but
%can probably do this based on regionprops. I'm going to do this on the
%background subtracted original green image of the HiK response.

%Also, ideally I'd like to erode the regions so that they are within the
%boundaries of the cell, to make sure that slight movements do not alter
%the mean intensity too much. 
D = bwdist(~L);
L2 = L;
L2(D < 6) = 0;
%figure, imshow(label2rgb(L, 'jet','k','shuffle'))
%figure, imshow(label2rgb(L2, 'jet','k','shuffle'))

statsL2 = regionprops(L2,imT,'Area','MeanIntensity','Centroid');
%statsL2 has a lot of NaN values where small labeled regions went to 0 in
%area, so there was no mean intensity to report.

%This part Andrew Paek helped me to design: remove cells above a certain
%max area. May have to adjust number based on histogram.

BigBadCells = find([statsL2.Area]>=7000);
for (j = 1:length(BigBadCells))
    L2(L2==BigBadCells(j))=0;
    [L3,num] = bwlabel(L2);
end
%figure, imshow(label2rgb(L3, 'jet','k','shuffle'))

%if there are NO big clumps, use this instead:
%[L3,num] = bwlabel(L2);
statsL3 = regionprops(L3,imT,'Area','MeanIntensity','Centroid');

%Now I will remove cells smaller than a certain area.
SmallBadCells = find([statsL3.Area]<=50);
for (x = 1:length(SmallBadCells));
    L3(L3==SmallBadCells(x))=0;
    [L4,num] = bwlabel(L3);
end
statsL4 = regionprops(L4,imT,'Area','MeanIntensity','Centroid');
figure, imshow(label2rgb(L4, 'jet','k','shuffle'))

%This part labels the L4 image with the numbers of each region
hold on
for k = 1:numel(statsL4)
    c = statsL4(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', 'white', 'fontsize',8);
end

%print(figure,[FileLoc '\' 'ROI.tif']); 