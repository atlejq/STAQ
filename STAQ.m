clear

config.basepath = 'C:\F\astro\matlab\m1test\';
config.darkPathRGB = 'parameters\darkframe10.tif';
config.darkPathH = 'parameters\darkframe20.tif';
config.filter = 'R';
config.align = 'R';

config.inputFormat = '.png';
config.maxStars = 10;
config.discardPercentile = 0.0;
config.medianOver = 40;
config.topMatchesMasterAlign = 5;
config.topMatchesMonoAlign = 5;

config.analyzeFrames = 0;  
config.findStackParameters = 1;
config.stackImages = 0;

ROI_y = 1:2822;
ROI_x = 1:4144;

e = 0.0005;

xvec={};
yvec={};
qual=[];

if(config.analyzeFrames == 1)      
    fileNameArray = getFileNames(config);
    
    for n = 1:height(fileNameArray)   
        lightFrame = imread(fileNameArray(n,:));
        lightFrame = lightFrame(ROI_y,ROI_x);
        starMatrix = analyzeStarField(lightFrame, config);
        stars(n) = length(starMatrix);
        background(n) = sum(sum(lightFrame));
        
        corrMatrix = flipud(sortrows(starMatrix',3));
        corrMatrix = corrMatrix';
        
        if(length(corrMatrix)>config.maxStars )
            corrMatrix = corrMatrix(1:2,1:config.maxStars);
        end
              
        xvec(n,:) = {corrMatrix(1,:)};
        yvec(n,:) = {corrMatrix(2,:)};
        
        if(mod(n,10)==0)
            n
        end
    end
    qual = stars/max(stars);
    
    [maxqual, q] = max(qual);
    refVectorX = xvec{q};
    refVectorY = yvec{q};
    maxQualFramePath = fileNameArray(q,:);
      
    save([config.basepath,'parameters\xvec',config.filter,'.mat'], 'xvec');
    save([config.basepath,'parameters\yvec',config.filter,'.mat'], 'yvec');
    save([config.basepath,'parameters\background',config.filter, '.mat'], 'background');  
    save([config.basepath,'parameters\qual',config.filter, '.mat'], 'qual');  
    save([config.basepath,'parameters\maxQualFramePath',config.filter, '.mat'], 'maxQualFramePath');
    save([config.basepath,'parameters\refVectorX',config.filter, '.mat'], 'refVectorX');
    save([config.basepath,'parameters\refVectorY',config.filter, '.mat'], 'refVectorY');
end

if(config.findStackParameters == 1)      
    xvec = importdata([config.basepath,'parameters\xvec',config.filter,'.mat']);
    yvec = importdata([config.basepath,'parameters\yvec',config.filter,'.mat']);
    background = importdata([config.basepath,'parameters\background',config.filter,'.mat']);
    qual = importdata([config.basepath,'parameters\qual',config.filter,'.mat']);
    maxQualFramePath = importdata([config.basepath,'parameters\maxQualFramePath',config.align,'.mat']);
    refVectorX = importdata([config.basepath,'parameters\refVectorX',config.filter,'.mat']);
    refVectorY = importdata([config.basepath,'parameters\refVectorY',config.filter,'.mat']);
    refVectorXAlign = importdata([config.basepath,'parameters\refVectorX',config.align,'.mat']);
    refVectorYAlign = importdata([config.basepath,'parameters\refVectorY',config.align,'.mat']);  
    
    qt = prctile(qual,config.discardPercentile*100);
    selectedFrames = find(qt<=qual);
    
    dx = zeros(1, length(selectedFrames));
    dy = zeros(1, length(selectedFrames));
    th = zeros(1, length(selectedFrames));
    discardFrames = uint32(zeros(1, length(selectedFrames)));  
    
    refTriangles = triangles(refVectorX,refVectorY);    
    refTrianglesAlign = triangles(refVectorXAlign,refVectorYAlign);

    refTriangles = sortrows(refTriangles,4);
    refTrianglesAlign = sortrows(refTrianglesAlign,4);

    [mth, mt, bla] = alignFrames(refVectorXAlign,refVectorYAlign,refTrianglesAlign,config.topMatchesMasterAlign,refVectorX,refVectorY,e);
    
    for i = 1:length(selectedFrames)           
       [theta, t, d] =  alignFrames(refVectorX,refVectorY,refTriangles,config.topMatchesMonoAlign,xvec{selectedFrames(i)},yvec{selectedFrames(i)},e);
       tmp = [cos(mth) -sin(mth); sin(mth) cos(mth)]*[t(1);t(2)] + [mt(1); mt(2)];    
       dx(i) = tmp(1);
       dy(i) = tmp(2);
       th(i) = theta+mth;
       discardFrames(i) = d;
    end 
   
    dx(discardFrames==1) = [];
    dy(discardFrames==1) = [];
    th(discardFrames==1) = [];
    selectedFrames(discardFrames==1) = [];
    
    maxQualFrame = imread(maxQualFramePath); 
    maxQualFrame = maxQualFrame(ROI_y,ROI_x);
    figure(1)
    imshow(maxQualFrame, 'Border', 'tight')
    hold on;
    plot(refVectorXAlign, refVectorYAlign, 'ro', 'MarkerSize',10); 
    for i=1:length(selectedFrames)
        R = [cos(th(i)) -sin(th(i)); sin(th(i)) cos(th(i))];
        t = [dx(i);dy(i)];
        debugMatrix = (R*[xvec{selectedFrames(i)}; yvec{selectedFrames(i)}]) + repmat(t, 1, length([xvec{selectedFrames(i)}; yvec{selectedFrames(i)}]));      
        plot(debugMatrix(1,:), debugMatrix(2,:), 'go');
        %plot(xvec{selectedFrames(i)},yvec{selectedFrames(i)}, 'go'); %Debug tracking errors
    end
    figure(2)
    plot(qual)
    hold on;
    plot(background/max(background))
    legend('Quality', 'Background')
    
    figure(3)
    plot(qual(selectedFrames))
    hold on;
    plot(background(selectedFrames)/max(background(selectedFrames)))
    legend('Quality', 'Background')
    
    save([config.basepath,'parameters\dx',config.filter,'.mat'],'dx');
    save([config.basepath,'parameters\dy',config.filter,'.mat'],'dy');
    save([config.basepath,'parameters\th',config.filter,'.mat'],'th');
    save([config.basepath,'parameters\selectedFrames',config.filter,'.mat'],'selectedFrames');
end

if(config.stackImages == 1)   
    load([config.basepath,'parameters\dx',config.filter,'.mat']);
    load([config.basepath,'parameters\dy',config.filter,'.mat']);    
    load([config.basepath,'parameters\th',config.filter,'.mat']);
    load([config.basepath,'parameters\selectedFrames',config.filter,'.mat'],'selectedFrames');
    load([config.basepath,'parameters\background',config.filter,'.mat'],'background');
    
    if(config.filter == "H")
        darkPath = config.darkPathH;
    else
        darkPath = config.darkPathRGB;
    end
        
    darkFrame = imread([config.basepath, darkPath]);
    darkFrame = rgb2gray(darkFrame);
    darkFrame = darkFrame(ROI_y,ROI_x);
    darkFrame = im2double(darkFrame);

    imarray = zeros(length(ROI_y),length(ROI_x));
    temparray = zeros(length(ROI_y),length(ROI_x));
    
    count = 1;
    tempcount = 1;
    fileNameArray = getFileNames(config);

    for i=1:length(selectedFrames)
        lightFrame = imread(fileNameArray(selectedFrames(i),:));
        lightFrame = lightFrame(ROI_y,ROI_x);
        lightFrame = im2double(lightFrame);      
        lightFrame = lightFrame*max(background(selectedFrames))/background(selectedFrames(i));
        lightFrame = lightFrame - darkFrame;
        
        trans = affine2d([cos(th(i)) sin(th(i)) 0; -sin(th(i)) cos(th(i)) 0; dx(i) , dy(i), 1]);
        outputView = imref2d([ceil(length(ROI_y)), ceil(length(ROI_x))]);
        lightFrame = imwarp(lightFrame,trans,'OutputView',outputView); 
        temparray(:,:,tempcount) = lightFrame(1:length(ROI_y),1:length(ROI_x));
        tempcount = tempcount + 1;
        if(mod(i,config.medianOver)==0)
            imarray(:,:,count) = median(temparray,3);
            clear temparray;
            count = count+1;
            tempcount = 1;
            i
        end
    end

    stackFrame = zeros(length(imarray(:,1,:)),length(imarray(1,:,:)));
    L = length(imarray(1,:,:));
    for i=1:length(imarray(:,1,:))
        for j=1:L
            testPixel = sort(imarray(i,j,:));
            stackFrame(i,j) = median(testPixel);
        end
    end
    
    imshow(stackFrame*5, 'Border', 'tight')    
    outputFrame = uint32(stackFrame*2^32);
    outputPath = [config.basepath, 'out\', num2str(length(selectedFrames)), '_', config.filter, '.tif'];
    save32BitImage(outputFrame, outputPath); 
end

function fileNameArray = getFileNames(config)
    fileNameArray = [];
    folders = dir([config.basepath, 'lights', '\','*']);
    for n = 3:size(folders,1)
        if(contains(folders(n).name, config.filter))
            lightPath = [config.basepath, 'lights', '\', folders(n).name, '\'];
            fileList = dir([lightPath, '*', config.inputFormat]);
                for m = 1:length(fileList)
                    path = [fileList(m).folder, '\', fileList(m).name];
                    fileNameArray = [fileNameArray; path];
                end
        end
    end       
end

function starMatrix = analyzeStarField(lightFrame, config)
    if(config.filter == "H")
        threshold = 0.5;
        factor = 3;
    else
        threshold = 0.9;
        factor = 1;
    end
    filteredImage = medfilt2(factor*lightFrame, [3 3]); 
    BW = imbinarize(filteredImage, threshold);
    [label_image] = bwlabel(BW,8);
    stats=regionprops(label_image,'BoundingBox'); 
    boundingBox = [stats.BoundingBox];
    starMatrix = [boundingBox(1:4:end)+boundingBox(3:4:end)/2; boundingBox(2:4:end)+boundingBox(4:4:end)/2;(boundingBox(3:4:end).^2+boundingBox(4:4:end).^2).^0.5];
end

function triangleParameters = triangles(x, y)
    triangleParameters = [];
    for i=1:length(x)-2
        for j=i+1:length(x)-1
            for k=j+1:length(x)            
                d = [sqrt((x(i)-x(j))^2 +(y(i)-y(j))^2), sqrt((x(j)-x(k))^2 +(y(j)-y(k))^2), sqrt((x(i)-x(k))^2 +(y(i)-y(k))^2)];          
                u = median(d)/max(d);
                v = min(d)/max(d);
                triangleParameters = [triangleParameters; [i j k u v]];            
            end
        end
    end
end

function [theta, t] = findRT(A, B)
    centroid_A = mean(A, 2);
    centroid_B = mean(B, 2);

    Am = A - repmat(centroid_A, 1, size(A, 2));
    Bm = B - repmat(centroid_B, 1, size(A, 2));

    H = Am * Bm';

    [U,S,V] = svd(H);
    R = V*U';

    if det(R) < 0
        V(:,2) = V(:,2) * -1;
        R = V*U';
    end
    
    theta = asin(R(2,1));
    t = -R*centroid_A + centroid_B;
end

function [theta, t, d] = alignFrames(refVectorX, refVectorY, refTriangles, topMatches, xvec, yvec, e)
        frameTriangles = triangles(xvec,yvec);
        vote = zeros(length(refVectorX),length(yvec), 'uint32');
        cVote = zeros(length(refVectorX),length(yvec), 'uint32'); 
          
        for a = 1:height(refTriangles)
            triangleList = find(((refTriangles(a,4)-e)<frameTriangles(:,4))&(frameTriangles(:,4)<(refTriangles(a,4)+e)));           
            for c = 1:length(triangleList)
               b=triangleList(c);               
                if (sqrt(((refTriangles(a,4) - frameTriangles(b,4))^2 + (refTriangles(a,5) - frameTriangles(b,5)))^2) < e)
                    vote((refTriangles(a,1)),(frameTriangles(b,1))) = vote(refTriangles(a,1),frameTriangles(b,1)) + 1;
                    vote((refTriangles(a,2)),(frameTriangles(b,2))) = vote(refTriangles(a,2),frameTriangles(b,2)) + 1;
                    vote((refTriangles(a,3)),(frameTriangles(b,3))) = vote(refTriangles(a,3),frameTriangles(b,3)) + 1;           
                end
            end   
        end
        for row = 1:height(vote)
            [maxRowVote, ind] = max(vote(row,:));
            cVote(row,ind) = maxRowVote-max(max(vote(row,[1:ind-1 ind+1:width(vote)])),max(vote([1:row-1 row+1:height(vote)],ind)));
        end
        
        cVote = max(cVote,0);
    
        [maxVote, maxVoteIndex] = max(cVote);
        
        votePairs = [1:1:length(maxVoteIndex); maxVoteIndex; maxVote];
        rankPairs = flipud(sortrows(votePairs',3));
        
        if(length(rankPairs(:,1))>=topMatches)           
            referenceMatrix = [refVectorX((rankPairs(1:topMatches,2)));refVectorY((rankPairs(1:topMatches,2)))];
            frameMatrix = [xvec(rankPairs(1:topMatches,1));yvec(rankPairs(1:topMatches,1))];                  
            [theta, t] = findRT(frameMatrix,referenceMatrix);
            d = 0;
        else
            disp(['Cannot stack frame']);
            theta = 0;
            t = 0;
            d = 1;
        end
end

function save32BitImage(outputFrame, outputPath)
    t = Tiff(outputPath,'w');
    tagstruct.ImageLength     = size(outputFrame,1);
    tagstruct.ImageWidth      = size(outputFrame,2);
    tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample   = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip    = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software        = 'MATLAB';
    t.setTag(tagstruct)
    t.write(outputFrame);
    t.close();
end