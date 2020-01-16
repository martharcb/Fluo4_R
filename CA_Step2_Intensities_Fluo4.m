%This is the second part of calcium imaging processing.
%Updated by Martha Bhattacharya on 11/09/19

%This section will take the labeled regions created during the first part
%(created on the brightest cells during the depolarization by High K+)
%and extract the mean intensity from those same regions across all of the
%images in a folder (one experimental trial).

%Note that you need to have run Step 1 immediately before this, and have
%the following in your workspace:
    %1) labeled image L4
    %2) numtiffs variable (number of tif images)
    %3) row dimensions of cleanstats (number of region labels)

prompt = {'What folder are the tif images in? Use full path:'};
dlgtitle = 'File Location';
dims = [1 50];
definput = {'O:\MBdata\Ca_imaging'};
FileLoc = inputdlg(prompt,dlgtitle,dims,definput);

%convert to a character output
FileLoc = char(FileLoc);

%make list of file names; create an empty struct for data for later
ImageDir = dir(FileLoc);
%exp_statsG = struct;
%exp_statsR = struct;
exp_statsF4 = zeros(numtiffs,size(statsL4,1));

NumFiles = size(ImageDir,1);
%ImageDir = ImageDir(1:NumFiles,:);

%For all the tif files, find the name, read in the image, background
%subtract with the rolling ball method used in Step 1.
for i = 1:NumFiles
    ImName = ImageDir(i).name;
    istif = regexp(ImName,'.tif');
    if (istif)>0
        disp(append('Analyzing frame ',num2str(i))) %reports frame
        fullname = [FileLoc '\' ImName];
        im = imread(fullname);
        imT = RollingBall(im,5,50);
        %imT = RollingBall(im(:,:,2),10,100);
%Here I am getting area and mean (all regions including tiny ones),
%removing NaN values and areas under 10 pixels, then just keeping the
%intensity (throwing out area) and transposing it.
%Do separately for green image, then red image.
       
          tempstatsF4 = regionprops(L4,imT,'MeanIntensity');
          %cleantempstatsF4 = rmmissing(struct2table(tempstatsF4));
          %cleantempstatsF4 = cleantempstatsF4(cleantempstatsF4.Area >= 10,:);
          %tempstatsF4 = tempstatsF4.MeanIntensity;
          tempstatsF4 = transpose(struct2array(tempstatsF4));
 %Add these data to the preallocated matrix in the proper row i (= frame)
            exp_statsF4(i,:) = tempstatsF4;
    end
end

writematrix(exp_statsF4, [FileLoc '\' 'stats.csv']);
%would like to use the file location variable to put this in the same
%directory as all of the tif images, eventually
