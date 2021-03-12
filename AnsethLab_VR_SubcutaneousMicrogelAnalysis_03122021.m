clear all; clc; tic; cd(userpath);
%% Editables %%

Folder = 'D:\OneDrive - UCB-O365\Work Stuff\Collaboration Files and Scripts\Anseth Lab - Varsha Rao Microgel Analysis\Widefield Images\4Channel\';
%Where are the images located? IMPORTANT: Use format 'FILEPATH\'. The apostrophes and ending slash are critical.

ROI_AssistedDraw = 1; %Do you want MatLab to help you contour your ROI along the underlying image? (1=yes,0=freehand)
    %In the Assisted mode, you can click to add "anchor points" to the ROI.
    %Hover the mouse near potential anchor points to allow MatLab to fit
    %the contour line to the underlying image. The contour-fit process
    %takes a little longer for larger images than smaller ones.
    
    %In the Freehand mode, you click and drag the cursor into whatever
    %shape you want. There are no anchor points and you can't release the
    %mouse button after you start.

BandDiameter = 200; %How "thick" (microns) do you want each erosion iteration to be?

ImageSave_Mask = 1; %Do you want to save an image of the original image (including the user-input mask)? (1=yes,0=no)
ImageSave_Erosion = 1; %Do you want to save an image showing the "erosion" outlines and cell area segmentation? (1=yes, 0=no)
ImageSave_ChannelMasks = 1; %Do you want to save images showing the cell segmentation for each channel? (1=yes, 0=no)

Channels = 4; %How many fluorescent channels are in the image?
CellChannels = [1 2 3]; %Which channels are cell-associated? This is used for segmentation of cell area. Format: [1 2 3 4], [1 3], etc.
Ch_DAPI = 1; %Which channel corresponds to DAPI/Hoechst signal?

MeanMultiplier = 1.5; %How stringent should the cell segmentation be, based on the MIP of all cell-associated channels?
%%This is based on the mean signal intensity of the entire MIP image.
%%Higher multipler = more stringent. Default = 1.

MultiThreshValue_Ch1 = 1; %How stringent do you want the segmentation of signal in Ch1? (Values must be 1-3, Default = 1)
MultiThreshValue_Ch2 = 1; %How stringent do you want the segmentation of signal in Ch2? (Values must be 1-3, Default = 1)
MultiThreshValue_Ch3 = 1; %How stringent do you want the segmentation of signal in Ch3? (Values must be 1-3, Default = 1)
MultiThreshValue_Ch4 = 1; %How stringent do you want the segmentation of signal in Ch4? (Values must be 1-3, Default = 1)

%% Analysis Pre-Reqs and Metadata %%

cd(Folder);
srcFiles = dir('*.nd2');

START = 1;
FINISH = 1; %length(srcFiles);

%% Analysis %%

for f = START:FINISH
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
    Metadata(f,1).Raw = I{1,2};
    Metadata(f,1).Commas = strfind(Metadata(f,1).Raw,',');
    Metadata(f,1).PixIdx(1,1) = strfind(Metadata(f,1).Raw,'Global dCalibration=')+20;
    Metadata(f,1).PixIdx(2,1) = find(Metadata(f,1).Commas>Metadata(f,1).PixIdx(1,1),1);
    Metadata(f,1).PixIdx(3,1) = Metadata(f,1).Commas(Metadata(f,1).PixIdx(2,1))-1;
    MicronsPerPixel(f,1) = str2num(extractBetween(string(Metadata(f,1).Raw),strfind(Metadata(f,1).Raw,'Global dCalibration=')+20,Metadata(f,1).PixIdx(3,1)));
    ResY = size(I{1,1}{1,1},1);
    ResX = size(I{1,1}{1,1},2);
    ZPlanes = (length(I{1,1})/Channels);
    
    %% Parsing Channels and Creating Max Intensity Projection Image %%
    disp('Parsing Channels and Creating Max Intensity Projection Image...');
    if f > START, clearvars Image; else end;
    for m = 1:ZPlanes
        ChannelPlanes.Ch1(m,1) = 1+(Channels*m-Channels);
        Image.Ch1(:,:,m) = I{1,1}{ChannelPlanes.Ch1(m,1),1};
        ChannelPlanes.Ch2(m,1) = 2+(Channels*m-Channels);
        Image.Ch2(:,:,m) = I{1,1}{ChannelPlanes.Ch2(m,1),1};
    if Channels>2
        ChannelPlanes.Ch3(m,1) = 3+(Channels*m-Channels);
        Image.Ch3(:,:,m) = I{1,1}{ChannelPlanes.Ch3(m,1),1}; 
    else
        Image.Ch3 = zeros(1); 
    end
    if Channels>3
        ChannelPlanes.Ch4(m,1) = 4+(Channels*m-Channels);
        Image.Ch4(:,:,m) = I{1,1}{ChannelPlanes.Ch4(m,1),1}; 
    else
        Image.Ch4 = zeros(1); 
    end
    end
    
    Image.Ch1MIP = max(Image.Ch1,[],3);
    Image.Ch2MIP = max(Image.Ch2,[],3);
    if Channels>2, Image.Ch3MIP = max(Image.Ch3,[],3); else end
    if Channels>3, Image.Ch4MIP = max(Image.Ch4,[],3); else end
    
    Image.CellMIP = uint16(zeros(ResY,ResX));
    if sum(CellChannels(:) == 1) == 1
        Image.CellMIP = max(Image.CellMIP,Image.Ch1MIP); else end
    if sum(CellChannels(:) == 2) == 1
        Image.CellMIP = max(Image.CellMIP,Image.Ch2MIP); else end
    if sum(CellChannels(:) == 3) == 1
        Image.CellMIP = max(Image.CellMIP,Image.Ch3MIP); else end
    if sum(CellChannels(:) == 4) == 1
        Image.CellMIP = max(Image.CellMIP,Image.Ch4MIP); else end
    
    Image.CellMIPMean = mean(Image.CellMIP,'all');
    Image.BWCellMIP = Image.CellMIP>(Image.CellMIPMean*MeanMultiplier);
    Image.BWCellMIP = bwareaopen(Image.BWCellMIP,10,8);
    Image.BWCellMIP = imerode(imdilate(Image.BWCellMIP,strel('disk',3)),strel('disk',3));
    
     %% Generating Figure and Collecting User Input %%
    disp('Generating Figure and Collecting User Input...');
    
     figure;
    FIG = imshow(imadjust(imcomplement(Image.CellMIP)));
    
    if ROI_AssistedDraw == 1
        ROI = drawassisted(FIG,'Color','r','Smoothing',5);
    else
        ROI = drawfreehand('Color','r');
    end
    
    Question_Continue = questdlg('Continue analysis with this mask?','Mask Confirmation','Continue','Try Again','Continue');
    
    if Ch_DAPI == 1, DAPI = Image.Ch1MIP;
    elseif Ch_DAPI == 2, DAPI = Image.Ch2MIP;
    elseif Ch_DAPI == 3, DAPI = Image.Ch3MIP;
    elseif Ch_DAPI == 4, DAPI = Image.Ch4MIP;
    else warning('Ch_DAPI values must be between 1 and 4 (inclusive).');
    end
    
    if strcmp(Question_Continue,'Continue') == 0
        FIG = imshow(imadjust(imcomplement(DAPI)));
        if ROI_AssistedDraw ==1, ROI = drawassisted(FIG,'Color','r','Smoothing',5);
        else ROI = drawfreehand('Color','r');
        end
        disp('Mask ROI confirmed...');
    else
        disp('Mask ROI confirmed...');
    end
    Image.BW1 = createMask(ROI);
    Image.BW1Perim = bwperim(Image.BW1);
    imshow(imoverlay(imadjust(Image.CellMIP),imdilate(Image.BW1Perim,strel('disk',5))));
    
    if ImageSave_Mask == 1
        disp('Saving Mask Figure Image...');
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
    
    close all;
    
    %% Erosion Masking %%
    disp('Creating erosion band masks...');
    if f > START, clearvars Erosion; else end
    Erosion.MicronsPerPixel = MicronsPerPixel(f,1);
    Erosion.ErosionBandWidth = BandDiameter;
    Erosion.PixelsPerBand = round(Erosion.ErosionBandWidth/Erosion.MicronsPerPixel);
    Erosion.strelrad = round(Erosion.PixelsPerBand);

    Erosion.Image(:,:,1) = Image.BW1;
    rp = regionprops(Erosion.Image(:,:,1),'MinorAxisLength');
    Erosion.MinorAxisLength = rp.MinorAxisLength;
    Erosion.ErosionDistance = Erosion.MinorAxisLength/2;
    Erosion.ErosionIterations = round(Erosion.ErosionDistance/Erosion.strelrad);

    Erosion.Perim(:,:,1) = Image.BW1Perim;
    Erosion.AllBandsMontage = Erosion.Perim(:,:,1);
    for e = 2:Erosion.ErosionIterations
        Erosion.Image(:,:,e) = imerode(Erosion.Image(:,:,e-1),strel('disk',Erosion.strelrad));
        Erosion.Perim(:,:,e) = bwperim(Erosion.Image(:,:,e));
        Erosion.AllBandsMontage = Erosion.AllBandsMontage + Erosion.Perim(:,:,e);
    end
    for n = 1:(Erosion.ErosionIterations-1)
        Erosion.BandMask(:,:,n) = Erosion.Image(:,:,n) - Erosion.Image(:,:,n+1);
        if n == Erosion.ErosionIterations-1
            Erosion.BandMask(:,:,n+1) = Erosion.Image(:,:,n+1);
        else
        end
    end
    Erosion.BandMask = logical(Erosion.BandMask);
    
    Erosion.MicrogelCellImage = Image.CellMIP.*uint16(Image.BW1);
    Erosion.MicrogelCellMask = and(Image.BWCellMIP,Image.BW1);
    
    Erosion.MicrogelCh1Image = Image.Ch1MIP.*uint16(Image.BW1);
    ChannelSeg.Ch1NZL = Erosion.MicrogelCh1Image(Erosion.MicrogelCh1Image~=0);
    ChannelSeg.Ch1MT = multithresh(ChannelSeg.Ch1NZL,3);
    Erosion.MicrogelCh1Mask = Erosion.MicrogelCh1Image>ChannelSeg.Ch1MT(MultiThreshValue_Ch1);
    Erosion.MicrogelCh1Mask = and(Erosion.MicrogelCh1Mask,Image.BWCellMIP);
    
    Erosion.MicrogelCh2Image = Image.Ch2MIP.*uint16(Image.BW1);
    ChannelSeg.Ch2NZL = Erosion.MicrogelCh2Image(Erosion.MicrogelCh2Image~=0);
    ChannelSeg.Ch2MT = multithresh(ChannelSeg.Ch2NZL,3);
    Erosion.MicrogelCh2Mask = Erosion.MicrogelCh2Image>ChannelSeg.Ch2MT(MultiThreshValue_Ch2);
    Erosion.MicrogelCh2Mask = and(Erosion.MicrogelCh2Mask,Image.BWCellMIP);
    
    if Channels > 2
       Erosion.MicrogelCh3Image = Image.Ch3MIP.*uint16(Image.BW1);
       ChannelSeg.Ch3NZL = Erosion.MicrogelCh3Image(Erosion.MicrogelCh3Image~=0);
       ChannelSeg.Ch3MT = multithresh(ChannelSeg.Ch3NZL,3);
       Erosion.MicrogelCh3Mask = Erosion.MicrogelCh3Image>ChannelSeg.Ch3MT(MultiThreshValue_Ch3);
       Erosion.MicrogelCh3Mask = and(Erosion.MicrogelCh3Mask,Image.BWCellMIP);
    else
    end
    
    if Channels > 3
       Erosion.MicrogelCh4Image = Image.Ch4MIP.*uint16(Image.BW1);
       ChannelSeg.Ch4NZL = Erosion.MicrogelCh4Image(Erosion.MicrogelCh4Image~=0);
       ChannelSeg.Ch4MT = multithresh(ChannelSeg.Ch4NZL,3);
       Erosion.MicrogelCh4Mask = Erosion.MicrogelCh4Image>ChannelSeg.Ch4MT(MultiThreshValue_Ch4);
       Erosion.MicrogelCh4Mask = and(Erosion.MicrogelCh4Mask,Image.BWCellMIP);
    else
    end
    
    Erosion.AllBandsMontage = imdilate(Erosion.AllBandsMontage,strel('disk',2));

    figimage1 = imoverlay(imadjust(Image.CellMIP),Erosion.AllBandsMontage);
    figimage2 = imoverlay(imadjust(Image.CellMIP),Erosion.MicrogelCellMask,[1 1 0]);

    C = [figimage1 figimage2];
    figure, imshow(C);

    if ImageSave_Erosion == 1
        disp('Saving Erosion Figure Image...');
        if f == START
            cd(Folder); cd Analysis; mkdir('ErosionImages'); cd ErosionImages;
        else cd(Folder); cd Analysis; cd ErosionImages;
        end
        ax = gca;
        slashfind = strfind(filename,'\');
        FigName_Erosion = append(filename((slashfind(end)+1):end-4),' Erosion Band Mask.jpg');

        exportgraphics(ax,FigName_Erosion,'ContentType','image','Resolution','400');
    else
    end
    
    
    %% Erosion Signal Analysis %%
    disp('Extracting signal information within each erosion band for each channel...');
    for m = 1:Erosion.ErosionIterations
        Erosion.CellAreaBandMask(:,:,m) = and(Erosion.BandMask(:,:,m),Erosion.MicrogelCellMask);
        ErosionProps(f,1).TotalCellArea(m,1) = sum(and(Erosion.BandMask(:,:,m),Erosion.MicrogelCellMask),'all');
        ErosionProps(f,1).BandArea(m,1) = sum(Erosion.BandMask(:,:,m),'all');
        ErosionProps(f,1).PercentCellAreaPerBand(m,1) = ErosionProps(f,1).TotalCellArea(m,1)/ErosionProps(f,1).BandArea(m,1)*100;
        
        ErosionImages(f,1).Ch1(:,:,m) = Image.Ch1MIP.*uint16(Erosion.BandMask(:,:,m));
        ErosionImages(f,1).Ch1CellImage(:,:,m) = Image.Ch1MIP.*uint16(Erosion.CellAreaBandMask(:,:,m));
        ErosionImages(f,1).Ch1CellMask(:,:,m) = and(Erosion.CellAreaBandMask(:,:,m),Erosion.MicrogelCh1Mask);
        ErosionProps(f,1).Ch1CellAreaSum(m,1) = sum(ErosionImages(f,1).Ch1CellMask(:,:,m),'all');
        ErosionProps(f,1).Ch1PercentofTotalCellArea(m,1) = ErosionProps(f,1).Ch1CellAreaSum(m,1)/ErosionProps(f,1).TotalCellArea(m,1)*100;
        ErosionProps(f,1).Ch1MeanBandIntensity(m,1) = mean(Image.Ch1MIP(Erosion.BandMask(:,:,m)),'all');
        ErosionProps(f,1).Ch1MeanCellIntensity(m,1) = mean(Image.Ch1MIP(Erosion.CellAreaBandMask(:,:,m)),'all');
        
        ErosionImages(f,1).Ch2(:,:,m) = Image.Ch2MIP.*uint16(Erosion.BandMask(:,:,m));
        ErosionImages(f,1).Ch2CellImage(:,:,m) = Image.Ch2MIP.*uint16(Erosion.CellAreaBandMask(:,:,m));
        ErosionImages(f,1).Ch2CellMask(:,:,m) = and(Erosion.CellAreaBandMask(:,:,m),Erosion.MicrogelCh2Mask);
        ErosionProps(f,1).Ch2CellAreaSum(m,1) = sum(ErosionImages(f,1).Ch2CellMask(:,:,m),'all');
        ErosionProps(f,1).Ch2PercentofTotalCellArea(m,1) = ErosionProps(f,1).Ch2CellAreaSum(m,1)/ErosionProps(f,1).TotalCellArea(m,1)*100;
        ErosionProps(f,1).Ch2MeanBandIntensity(m,1) = mean(Image.Ch2MIP(Erosion.BandMask(:,:,m)),'all');
        ErosionProps(f,1).Ch2MeanCellIntensity(m,1) = mean(Image.Ch2MIP(Erosion.CellAreaBandMask(:,:,m)),'all');
        
        if Channels>2
            ErosionImages(f,1).Ch3(:,:,m) = Image.Ch3MIP.*uint16(Erosion.BandMask(:,:,m));
            ErosionImages(f,1).Ch3CellImage(:,:,m) = Image.Ch3MIP.*uint16(Erosion.CellAreaBandMask(:,:,m));
            ErosionImages(f,1).Ch3CellMask(:,:,m) = and(Erosion.CellAreaBandMask(:,:,m),Erosion.MicrogelCh3Mask);
            ErosionProps(f,1).Ch3CellAreaSum(m,1) = sum(ErosionImages(f,1).Ch3CellMask(:,:,m),'all');
            ErosionProps(f,1).Ch3PercentofTotalCellArea(m,1) = ErosionProps(f,1).Ch3CellAreaSum(m,1)/ErosionProps(f,1).TotalCellArea(m,1)*100;
            ErosionProps(f,1).Ch3MeanBandIntensity(m,1) = mean(Image.Ch3MIP(Erosion.BandMask(:,:,m)),'all');
            ErosionProps(f,1).Ch3MeanCellIntensity(m,1) = mean(Image.Ch3MIP(Erosion.CellAreaBandMask(:,:,m)),'all');
        else
        end
        
        if Channels>3
            ErosionImages(f,1).Ch4(:,:,m) = Image.Ch4MIP.*uint16(Erosion.BandMask(:,:,m));
            ErosionImages(f,1).Ch4CellImage(:,:,m) = Image.Ch4MIP.*uint16(Erosion.CellAreaBandMask(:,:,m));
            ErosionImages(f,1).Ch4CellMask(:,:,m) = and(Erosion.CellAreaBandMask(:,:,m),Erosion.MicrogelCh4Mask);
            ErosionProps(f,1).Ch4CellAreaSum(m,1) = sum(ErosionImages(f,1).Ch4CellMask(:,:,m),'all');
            ErosionProps(f,1).Ch4PercentofTotalCellArea(m,1) = ErosionProps(f,1).Ch4CellAreaSum(m,1)/ErosionProps(f,1).TotalCellArea(m,1)*100;
            ErosionProps(f,1).Ch4MeanBandIntensity(m,1) = mean(Image.Ch4MIP(Erosion.BandMask(:,:,m)),'all');
            ErosionProps(f,1).Ch4MeanCellIntensity(m,1) = mean(Image.Ch4MIP(Erosion.CellAreaBandMask(:,:,m)),'all');
        else
        end
    end
    
    if ImageSave_ChannelMasks == 1
        disp('Generating and saving Channel Mask Images...');
        figimagech1 = imoverlay(imadjust(Image.Ch1MIP),Erosion.MicrogelCh1Mask);
        figimagech2 = imoverlay(imadjust(Image.Ch2MIP),Erosion.MicrogelCh2Mask);
        if Channels > 2
            figimagech3 = imoverlay(imadjust(Image.Ch3MIP),Erosion.MicrogelCh3Mask);
        else
        end
        if Channels > 3
            figimagech4 = imoverlay(imadjust(Image.Ch4MIP),Erosion.MicrogelCh4Mask);
        else
        end
        
        if Channels == 2
            D = [figimagech1 figimagech2];
        elseif Channels == 3
            D = [figimagech1 figimagech2 figimagech3];
        elseif Channels == 4
            D = [figimagech1 figimagech2 figimagech3 figimagech4];
        else
        end
        
        figure,imshow(D); title('Channel Segmentation Masks');

        if f == START
            cd(Folder); cd Analysis; mkdir('ChannelMaskImages'); cd ChannelMaskImages;
        else cd(Folder); cd Analysis; cd ChannelMaskImages;
        end
        ax = gca;
        slashfind = strfind(filename,'\');
        FigName_ChannelSeg = append(filename((slashfind(end)+1):end-4),' Channel Segmentation Mask.jpg');
        exportgraphics(ax,FigName_ChannelSeg,'ContentType','image','Resolution','400');
    else
    end   
    
    close all;
end

 disp('Collating results and saving to .mat file...');
 cd(Folder); cd Analysis; mkdir('ErosionAnalysisResults'); cd ErosionAnalysisResults;
 save('ErosionAnalysisResults.mat','ErosionProps','-v7.3');
    