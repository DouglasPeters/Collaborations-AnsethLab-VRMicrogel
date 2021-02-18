clear all; clc; tic; cd(userpath);
%% Editables %%

Folder = 'C:\Users\dplad\OneDrive - UCB-O365\Work Stuff\Collaboration Files and Scripts\Anseth Lab - Varsha Rao Microgel Analysis\';
%Where are the images located? IMPORTANT: Use format 'FILEPATH\'. The apostrophes and ending slash are critical.

ROI_Assisted = 1; %Do you want MatLab to help you contour your ROI along the underlying image? (1=yes,0=freehand)
Analysis_Erosion = 1; %Do you want to perform Erosion analysis on the mask? (1=yes,0=no)
    Analysis_ErosionBandWidth = 500; %How "thick" (microns) do you want each erosion iteration to be?

ImageSave_Mask = 0; %Do you want to save an image of the original image (including the user-input mask)? (1=yes,0=no)
ImageSave_Erosion = 0; %Do you want to save an image showing the "erosion" outlines? (1=yes, 0=no)

Channels = 2; %How many fluorescent channels are in the image?

CH_DAPI = 1; %Which channel corresponds to DAPI signal?

%% Analysis Pre-Reqs and Metadata %%

cd(Folder);
srcFiles = dir('*.nd2');

START = 1;
FINISH = length(srcFiles);

figure;

%% Analysis %%

for f = START:1 %FINISH
    time(f,1).ElapsedSeconds = toc;
    
    clc
    filename = strcat(Folder,srcFiles(f).name);
    progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2);
    
    if f == START
        cd(Folder); mkdir('Analysis'); cd(Folder);
    else end
    
    if progress < 10
        disp('Estimated time remaining will display after 10% of images are analyzed...');
    else
        time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/((FINISH-START+1)-(FINISH-f));
        time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH-START+1);
        time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
        time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
        time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
        estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
        estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
        disp(estimate);
        disp(estimate2);
    end

    Results(f).FileName = srcFiles(f).name;
    
    I = bfopen(filename);
    Metadata = I{1,2};
    MicronsPerPixel = str2num(extractBetween(string(Metadata),strfind(Metadata,'Global dCalibration=')+20,strfind(Metadata,'Global dCalibration=')+26));
    ResY = size(I{1,1}{1,1},1);
    ResX = size(I{1,1}{1,1},2);
    ZPlanes = (length(I{1,1})/Channels);
    
    %% Parsing Channels and Creating Max Intensity Projection Image %%
    disp('Parsing Channels and Creating Max Intensity Projection Image...');
    for m = 1:ZPlanes
        Image.Ch1(:,:,m) = I{1,1}{m,1}(:,:);
        Image.Ch2(:,:,m) = I{1,1}{ZPlanes+m,1}(:,:);
    if Channels>2, Image.Ch3(:,:,m) = I{1,1}{(2*ZPlanes)+m,1}(:,:); 
    else Image.Ch3 = zeros(1); end
    if Channels>3, Image.Ch4(:,:,m) = I{1,1}{(3*ZPlanes)+m,1}(:,:); 
    else Image.Ch4 = zeros(1); end
    end
    
    Image.Ch1MIP = max(Image.Ch1(:,:),3);
    Image.Ch2MIP = max(Image.Ch2(:,:),3);
    if Channels>2, Image.Ch3MIP = max(Image.Ch3(:,:),3); else end
    if Channels>3, Image.Ch4MIP = max(Image.Ch4(:,:),3); else end
    
    %% Generating Figure and Collecting User Input %%
    disp('Generating Figure and Collecting User Input...');
    
    if CH_DAPI == 1, DAPI = Image.Ch1MIP;
    elseif CH_DAPI == 2, DAPI = Image.Ch2MIP;
    elseif CH_DAPI == 3, DAPI = Image.Ch3MIP;
    elseif CH_DAPI == 4, DAPI = Image.Ch4MIP;
    else
    end
    
    DAPI_MeanIntensity = mean(DAPI,'all');
    
    FIG = imshow(imadjust(imcomplement(DAPI)));
    
    if ROI_Assisted == 1
        ROI = drawassisted(FIG,'Color','r','Smoothing',5);
    else
        ROI = drawfreehand('Color','r');
    end
    
    Question_Continue = questdlg('Continue analysis with this mask?','Mask Confirmation','Continue','Try Again','Continue');
    
    if strcmp(Question_Continue,'Continue') == 0 
        FIG = imshow(imadjust(imcomplement(DAPI)));
        if ROI_Assisted ==1, ROI = drawassisted(FIG,'Color','r','Smoothing',5);
        else ROI = drawfreehand('Color','r');
        end
        disp('Mask ROI confirmed...');
    else
        disp('Mask ROI confirmed...');
    end
    Image.BW1 = createMask(ROI);
    Image.BW1Perim = bwperim(Image.BW1);
    imshow(imoverlay(imadjust(DAPI),imdilate(Image.BW1Perim,strel('disk',5))));
    
    if ImageSave_Mask == 1
        disp('Saving Figure Image...');
        if f == START
            cd(Folder); cd Analysis; mkdir('MaskImages'); cd MaskImages;
        else cd(Folder); cd Analysis; cd MaskImages;
        end
        ax = gca;
        slashfind = strfind(filename,'\');
        FigName_Mask = append(filename((slashfind(end)+1):end-4),' ROI Mask.jpg');
        
        exportgraphics(ax,FigName_Mask,'ContentType','image','Resolution','400');
    else
    end
    
    %% Erosion Analysis %%
if Analysis_Erosion == 1
    disp('Performing erosion analysis...');
    Erosion.MicronsPerPixel = MicronsPerPixel;
    Erosion.ErosionBandWidth = Analysis_ErosionBandWidth;
    Erosion.PixelsPerBand = round(Erosion.ErosionBandWidth/Erosion.MicronsPerPixel);
    Erosion.strelrad = round(Erosion.PixelsPerBand/2);

    Erosion.Image(:,:,1) = Image.BW1;
    rp = regionprops(Erosion.Image(:,:,1),'MinorAxisLength');
    Erosion.MinorAxisLength = rp.MinorAxisLength;
    Erosion.ErosionDistance = Erosion.MinorAxisLength/2;
    Erosion.ErosionIterations = round(Erosion.ErosionDistance/Erosion.strelrad);
  
    Erosion.Perim(:,:,1) = Image.BW1Perim;
    Erosion.Montage = Erosion.Perim(:,:,1);
    for e = 2:Erosion.ErosionIterations
        Erosion.Image(:,:,e) = imerode(Erosion.Image(:,:,e-1),strel('disk',Erosion.strelrad));
        Erosion.Perim(:,:,e) = bwperim(Erosion.Image(:,:,e));
        Erosion.Montage = Erosion.Montage + Erosion.Perim(:,:,e);
    end
    for n = 1:(Erosion.ErosionIterations-1)
        Erosion.Band(:,:,n) = Erosion.Image(:,:,n) - Erosion.Image(:,:,n+1);
        if n == Erosion.ErosionIterations-1
            Erosion.Band(:,:,n+1) = Erosion.Image(:,:,n+1);
        else
        end
    end
    Erosion.Band = logical(Erosion.Band);
    
    Erosion.Microgel = DAPI.*uint16(Image.BW1);
    Erosion.CellAreaMask = Erosion.Microgel>DAPI_MeanIntensity;
    Erosion.Montage = imdilate(Erosion.Montage,strel('disk',2));
    
    figimage1 = imoverlay(imadjust(DAPI),Erosion.Montage);
    figimage2 = imoverlay(imadjust(DAPI),Erosion.CellAreaMask);
    
    C = [figimage1 figimage2];
    imshow(C);

    for m = 1:Erosion.ErosionIterations
        ErosionProps(f,1).CellArea(m,1) = sum(and(Erosion.Band(:,:,m),Erosion.CellAreaMask),'all');
        ErosionProps(f,1).BandArea(m,1) = sum(Erosion.Band(:,:,m),'all');
        ErosionProps(f,1).PercentCellAreaPerBand(m,1) = ErosionProps(f,1).CellArea(m,1)/ErosionProps(f,1).BandArea(m,1)*100;
        
        ErosionProps(f,1).Ch1Bands(:,:,m) = Image.Ch1MIP.*uint16(Erosion.Band(:,:,m));
        ErosionProps(f,1).Ch1MeanIntensity(m,1) = mean(Image.Ch1MIP(Erosion.Band(:,:,m)),'all');
        
        ErosionProps(f,1).Ch2Bands(:,:,m) = Image.Ch2MIP.*uint16(Erosion.Band(:,:,m));
        ErosionProps(f,1).Ch2MeanIntensity(m,1) = mean(Image.Ch2MIP(Erosion.Band(:,:,m)),'all');
        
        if Channels>2
            ErosionProps(f,1).Ch3Bands(:,:,m) = Image.Ch3MIP.*uint16(Erosion.Band(:,:,m));
            ErosionProps(f,1).Ch3MeanIntensity(m,1) = mean(Image.Ch3MIP(Erosion.Band(:,:,m)),'all');
        else
        end
        if Channels>3
            ErosionProps(f,1).Ch4Bands(:,:,m) = Image.Ch4MIP.*uint16(Erosion.Band(:,:,m));
            ErosionProps(f,1).Ch4MeanIntensity(m,1) = mean(Image.Ch4MIP(Erosion.Band(:,:,m)),'all');
        else
        end
    end
else
end
    
end
