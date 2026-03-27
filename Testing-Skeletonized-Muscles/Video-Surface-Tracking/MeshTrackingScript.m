
% % im = imread("Example Gradient.png");
% % im = rgb2gray(im);
% im2 = val_im;
% for i = 1:5
% im2 = imtophat(im2,strel('disk',5));
% end
% bw1 = im2 > 0.05;
% bw2 = bwareaopen(bw1,40);
% 
% figure()
% imshow(im)
% figure()
% imshow(val_im)
% figure()
% imshow(im2)
% figure()
% imshow(bw1)
% figure()
% imshow(bw2)

% im_og = im;
% for i = 1:10
%     bg = imopen(im,strel('disk',100));
%     im = im - bg;
% end
% figure()
% imshow(im_og)
% figure()
% imshow(bg)
% figure()
% imshow(im)

% vid_obj = VideoReader("2 Attachment Point Bar Inflation Videos\White Mesh\80mm Bar.MOV");
% 
% frame = read(vid_obj,1080);
% 
% % Isolate Mesh
% hsv_frame = rgb2hsv(frame);
% sat_mask = hsv_frame(:,:,2) < 0.1;
% val_im = hsv_frame(:,:,3);
% val_im(~sat_mask) = 0;
% val_mask = imregionalmax(val_im);
% 
% figure()
% imshow(frame)
% figure()
% imshow(sat_mask);
% figure()
% imshow(val_im)
% figure()
% imshow(val_mask)


vid_obj = VideoReader("2 Attachment Point Bar Inflation Videos\White Mesh\80mm Bar.MOV");

frame = read(vid_obj,1080);
binary = ~binarize_frame(frame,'method',"hsv match",'range',[1,0.12,0.35]);
binary = imfill(binary,'holes');
cleaned_binary1 = bw_clean_im(binary,10000);
cleaned_binary2 = bwmorph(cleaned_binary1,'thin',25);
cleaned_binary3 = bwmorph(cleaned_binary2,'majority');

dist = bwdist(cleaned_binary3);
dist_disp = dist/max(dist,[],'all');
ridge = imregionalmax(dist);
cleaned_ridge1 = imdilate(ridge,strel("square",50));
cleaned_ridge2 = bwareaopen(cleaned_ridge1, 75^2);
cleaned_ridge3 = bwskel(cleaned_ridge2);

% skel = bwskel(~cleaned_binary,'MinBranchLength',5000);
% 
% % convert skeleton to graph and find longest continuous path to treat as
% % midline
% endpts = bwmorph(skel,'endpoints');
% [y,x] = find(endpts);



figure()
imshow(frame);
figure()
imshow(binary);
figure()
imshow(cleaned_binary1);
figure()
imshow(cleaned_binary2);
figure()
imshow(cleaned_binary3);
figure()
imshow(dist_disp);
figure()
imshow(ridge);
figure()
imshow(cleaned_ridge1);
figure()
imshow(cleaned_ridge2);
figure()
imshow(cleaned_ridge3);
figure()
imshow(labeloverlay(frame,imdilate(cleaned_ridge3,strel("disk",10,4)),'Transparency',0))
figure()
hold on
imshow(frame);
if ~isempty(points1)
    scatter(points1(:,1), points1(:,2), 'filled')
end
if ~isempty(points2)
    scatter(points2(:,1), points2(:,2), 'filled')
end
hold off

function cleaned_im = bw_clean_im(binary_frame, hole_size)
% Purpose:
%       - Takes a binary image of a McKibben muscle and cleans it by
%         removing small holes and stray regions of black
%
% Arguments:
%       - binary_frame: a binarized image of the McKibben muscle where the
%       mesh is black and the balloon and background are white
%       - hole_size : holes smaller than this will be closed
%
% Returns:
%       - cleaned_im: the cleaned binary image

    % check arguments
    arguments
        binary_frame (:,:) {mustBeLogical}
        hole_size (1,1) {mustBeUint8orDouble} = 50
    end

    cleaned_im = ~bwmorph(~binary_frame,'close');       % connect small gaps   
    cleaned_im = ~bwareaopen(~cleaned_im, hole_size);   % fill small holes
    cleaned_im = bwareaopen(cleaned_im, hole_size);     % fill small holes


end

function [mesh_crossings, gap_centers] = get_points(binary_frame, options)
% Purpose:
%       - Takes a binary image of a McKibben muscle and identifies the
%         positions of mesh strand intersections and/or the positions of
%         the centroids of the diamond-like gaps between mesh strands
%
% Arguments:
%       - binary_frame: a binarized image of the McKibben muscle where the
%       mesh is black and the balloon and background are white
%       - options:
%           - method: method used to isolate mesh
%               - "intersections": search for mesh crossings
%               - "holes": search for holes in mesh
%               - "both": search for both types of points
%           - hole_range: range of hole area in pixels allowed
%               - Holes smaller than this will be closed
%               - Holes larger than this will be ignored (for background)
%           - cluster_dist: points in a cluste closer than this together
%             will be grouped together to form one point in output
% Returns:
%       - mesh_crossings: a [Nx2] array of points at the mesh intersections
%       - gap_centers: a [Nx2] array of points at the centroid of mesh gaps

    % check arguments
    arguments
        binary_frame (:,:) {mustBeLogical}

        options.method (1,1) string ...
            {mustBeMember(options.method, ["intersections","holes","both"])} ...
            = "intersections"

        options.hole_range {mustBeRange, mustBeUint8orDouble} = [50,10000]
        options.cluster_dist (1,1) {mustBeScalar} = 10
    end

    % run image processing
    switch options.method
        case "intersections"
            mesh_crossings = get_intersections();
            gap_centers = [];
        
        case "holes"
            mesh_crossings = [];
            gap_centers = get_holes();

        case "both"
            mesh_crossings = get_intersections();
            gap_centers = get_holes();
    end

    function intersections = get_intersections()
        % get all points where skeleton branches
        skeleton = bwskel(~binary_frame,'MinBranchLength',100);
        junctions = bwmorph(skeleton,'branchpoints');
        [y,x] = find(junctions);
        pts = [x,y];

        if options.cluster_dist ~= 0
            intersections = cluster_points(pts);
        else
            intersections = pts;
        end
    end

    function holes = get_holes()
        connected_regions = bwconncomp(binary_frame);
        stats = regionprops(connected_regions, 'Centroid','Area');
        areas = [stats.Area];
        valid_points = areas >= options.hole_range(1) & areas <= options.hole_range(2);
        holes = vertcat(stats(valid_points).Centroid);
        
        if options.cluster_dist ~= 0
            holes = cluster_points(holes);
        end
    end

    function grouped_points = cluster_points(points)
        dists = pdist2(points,points);                      % get distances between every pair of points
        A = (dists < options.cluster_dist) & (dists > 0);   % find all distances less than cluster_dist
        G = graph(A);                                       % build a graph where points are connected by edges if they are less than 30 away from each other
        bins = conncomp(G);                                 % assign each point to a group with all the points its nearby to
        numClusters = max(bins);                            % number of distinct groups is number of true intersections
        grouped_points = zeros(numClusters,2);              % fill in center of clusters with average of each point per cluster
        for k = 1:numClusters
            grouped_points(k,:) = mean(points(bins == k,:), 1);
        end
    end
end

function bw_frame = binarize_frame(frame, options)
% Purpose:
%       - Takes a frame of video of McKibben contraction and binarizes it
%         to a black and white image where only the mesh is white
%
% Arguments:
%       - frame: the frame of video in RGB or greyscale (uint8 or float accepted)
%       - options:
%           - method: method used to isolate mesh
%               - "color match": uses RGB color to isolate mesh
%                   - use 0-255 domain for target and range
%               - "val match": does not use range, isolates mesh as all
%                  pixels above target value,
%                   - use 0-1 domain for target and range
%               - "hsv match: uses hue, sat, and value to isolate mesh,
%                   - use 0-1 domain for target and range
%           - target: color or val of the mesh in image
%           - range: +\- range value for isolation method chosen
%             (larger values include more of the image)
% Returns:
%       - bw_frame: binarized image where ideally only mesh is pure white

    % check arguments
    arguments
        frame {mustBeImage,mustBeUint8orDouble}

        options.method (1,1) string ...
            {mustBeMember(options.method, ["color match","val match","hsv match"])} ...
            = "color match"

        options.target = []
        options.range = []
    end

    % convert input image to double
    if isa(frame,'uint8')
        frame = im2double(frame);
    end

    % default values for different methods
    default_color = [241,20,20];
    default_val_threshold = 0.85;
    default_hsv = [1,0,1];
    default_color_range = 75;
    default_hsv_range = [0.3,0.04,0.5];

    % more complex argument handling for each method with binarization
    % performed at the end of each case via subfunctions
    switch options.method
        case "color match"
            mustBeRGB(frame)

            if isempty(options.target)
                options.target = default_color;
            end
            mustBeUint8orDouble(options.target);
            mustBe1x3(options.target);
            options.target = options.target/255;
        
            if isempty(options.range)
                options.range = default_color_range;
            end
            mustBeUint8orDouble(options.range);
            mustBeScalar(options.range);
            options.range = options.range/255;
            
            bw_frame = binarize_by_color();

        case "val match"
            if isempty(options.target)
                options.target = default_val_threshold;
            end
            mustBeUint8orDouble(options.target);
            mustBeScalar(options.target);

            bw_frame = binarize_by_val();
        
        case "hsv match"
            mustBeRGB(frame)

            if isempty(options.target)
                options.target = default_hsv;
            end
            mustBeUint8orDouble(options.target);
            mustBe1x3(options.target);

        
            if isempty(options.range)
                options.range = default_hsv_range;
            end
            mustBeUint8orDouble(options.range);
            mustBe1x3(options.range)

            bw_frame = binarize_by_hsv();
    end

    function mask = binarize_by_color()
        frame = double(frame);
        color = double(reshape(options.target,[1,1,3]));
        dist = sqrt(sum((frame - color).^2, 3));
        mask = dist < options.range;
    end

    function mask = binarize_by_val()
        if length(size(frame)) == 3
            frame = rgb2hsv(frame);
            mask = frame(:,:,3) > options.target;
        else
            mask = frame > options.target;
        end
    end

    function mask = binarize_by_hsv()
        frame = rgb2hsv(frame);

        % shift hue so that the target is 0.5
        shift_val = mod(0.5 - options.target(1),1);
        frame(:,:,1) = frame(:,:,1) + shift_val;
        frame(frame > 1) = frame(frame > 1) - 1;

        % binarize image by ranges
        hue_bounds = [0.5 - options.range(1),...
                      0.5 + options.range(1)];
        sat_bounds = [options.target(2) - options.range(2),...
                      options.target(2) + options.range(2)];
        val_bounds = [options.target(3) - options.range(3),...
                      options.target(3) + options.range(3)];
        mask = isbetween(frame(:,:,1), hue_bounds(1), hue_bounds(2)) & ...
                isbetween(frame(:,:,2), sat_bounds(1), sat_bounds(2)) & ...
                isbetween(frame(:,:,3), val_bounds(1), val_bounds(2));
    end
end

function mustBeRGB(x)
    dimensions = size(x);
    if length(dimensions) ~= 3
        error('Image must be RGB')
    elseif dimensions(3) ~= 3
        error('Color image must have 3 color channels')
    end
end

function mustBeImage(x)
    dimensions = size(x);
    if length(dimensions) == 3
        if dimensions(3) ~= 3
            error('Color image must have 3 color channels')
        end
    elseif length(dimensions) ~= 2
        error('Image must be either greyscale or color')
    end
end

function mustBeScalar(x)
    if ~isscalar(x)
        error('Input must be scalar')
    end
end

function mustBeUint8orDouble(x)
    if ~(isa(x,"uint8") || isa(x,"double"))
        error('Input must be uint8 or double type')
    end
end

function mustBe1x3(x)
    if size(x) ~= 3
        error('Input must be a 1x3 array')
    end
end

function mustBeLogical(x)
    if ~islogical(x)
        error('Input must be logical')
    end
end

function mustBeRange(x)
    if size(x) ~= 2
        error('Input must be a 1x2 array')
    end
end



% mesh_hue = 356/360;
% mesh_hue_range = 3/360;
% hue_bound = [mesh_hue - mesh_hue_range, mesh_hue + mesh_hue_range];
% 
% mesh_sat = 75/100;
% mesh_sat_range = 20/100;
% sat_bound = [mesh_sat - mesh_sat_range, mesh_sat + mesh_sat_range];
% 
% mesh_val = 95/100;
% mesh_val_range = 10/100;
% val_bound = [mesh_val - mesh_val_range, mesh_val + mesh_val_range];
% 
% reverseStr = '';        % Initialize a variable to hold the backspaces
% 
% color = [241,20,20];
% color = double(reshape(color,1,1,3));

% frame = readFrame(vid_obj);
% frame = double(frame);
% frame = frame(:,:,2);
% k0 = .7;
% k1 = k0*(100);
% k1 = ones(size(frame))*k1;
% frame = (atan(k0.*frame-k1) + pi/2) * 255/pi;
% frame = (frame.^k)/(255.^(k-1));
% frame = uint8(frame);

% frame = double(frame);
% dist = sqrt(sum((frame - color).^2, 3));
% mask = dist<140;
% frame(:,:,1) = mask.*frame(:,:,1);
% frame(:,:,2) = mask.*frame(:,:,2);
% frame(:,:,3) = mask.*frame(:,:,3);
% frame = uint8(frame);
% frame = frame(:,:,1);

% frame = frame(:,:,[3 2 1]); % swap r and b so that we're not looking for value of 1
% frame = rgb2hsv(frame);
% hue = frame(:,:,1);
% % hue = mod(hue + 0.5,1); % wrap hue around cuz 1 and 0 are almost the same
% % hue(hue == 0) = 1;      % number and full red (color of mesh) is 1
% sat = frame(:,:,2);
% val = frame(:,:,3);
% mask = val > 0.9;       % build mask out of value

%%%%

% frame = imread("Video Tracking Test Image.png");
% green = frame(:,:,2);
% binary = imbinarize(green);
% binary = imrotate(~binary,26);
% % binary = bwmorph(binary,'spur');
% % binary = bwmorph(binary,'thin',Inf);
% skeleton = bwskel(binary,'MinBranchLength',100);
% junctions = bwmorph(skeleton,'branchpoints');
% [y,x] = find(junctions);
% % v average close-by points to get true centers v
% pts = [x,y];
% dists = pdist2(pts,pts);            % get distances between every pair of points
% A = (dists < 30) & (dists > 0);     % find all distances less than 30
% G = graph(A);                       % build a graph where points are connected by edges if they are less than 30 away from each other
% bins = conncomp(G);                 % assign each point to a group with all the points its connected to
% numClusters = max(bins);            % number of distinct groups is number of true intersections
% centers = zeros(numClusters,2);     % fill in center of clusters with average of each point per cluster
% for k = 1:numClusters
%     centers(k,:) = mean(pts(bins == k,:), 1);
% end
% figure()
% imshow(frame);
% figure()
% imshow(green);
% figure()
% imshow(binary);
% figure()
% imshow(skeleton);
% figure(); hold on;
% imshow(binary);
% scatter(centers(:,1),centers(:,2),'filled');
% hold off;
% 
% 
% frame = imread("Video Tracking Test Image.png");
% green = frame(:,:,2);
% binary = imbinarize(green);
% binary = imrotate(binary,26);
% binary = bwmorph(binary,'fill',Inf);                % make sure no holes in balloon regions
% binary = bwareaopen(binary,1000);                   % delete regions smaller than 1000 pixels in area
% connected_regions = bwconncomp(binary);
% stats = regionprops(connected_regions, 'Centroid');
% centroids = vertcat(stats.Centroid);
% figure()
% imshow(frame)
% figure()
% imshow(green)
% figure()
% imshow(binary)
% figure(); hold on;
% imshow(binary)
% scatter(centroids(:,1), centroids(:,2), 'filled')


%%%%

% corners = detectHarrisFeatures(binary, 'MinQuality', .15);

% edges = edge(binary,'canny');

% % Calculate the Hough transform
% [H, theta, rho] = hough(edges);
% 
% % Find peaks in the Hough transform to get the most prominent lines
% % The second argument specifies the number of lines to detect (e.g., 5)
% peaks = houghpeaks(H, 100, 'threshold', ceil(0.1*max(H(:)))); 
% 
% % Extract the actual line segments from the edge image
% lines = houghlines(edges, theta, rho, peaks, 'FillGap', 25, 'MinLength', 30);



% figure()
% imshow(hue)
% figure()
% imshow(sat)
% figure()
% imshow(val)
% figure()
% imshow(mask);
% figure()
% imshow(edges);

% % Display the original image and hold the plot for drawing lines
% figure, imshow(frame), hold on;
% 
% % Iterate through the lines and plot them
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'red');
% 
%    % Plot beginnings and ends of lines
%    plot(xy(1,1), xy(1,2), 's', 'LineWidth', 2, 'Color', 'yellow');
%    plot(xy(2,1), xy(2,2), 's', 'LineWidth', 2, 'Color', 'green');
% end
% title('Lines found between edges');
% hold off;

%impixelinfo;

%frame = im2uint8(frame);

% frame = rgb2hsv(frame);
% mask = isbetween(frame(:,:,1), hue_bound(1), hue_bound(2)) & ...
%          isbetween(frame(:,:,2), sat_bound(1), sat_bound(2)) & ...
%          isbetween(frame(:,:,3), val_bound(1), val_bound(2));
% frame(:,:,3) = mask.*frame(:,:,3);
% frame = hsv2rgb(frame);
% frame(:,:,1) = frame(:,:,1).*mask;






%%
%{
vid_out_obj = VideoWriter('edited_vid.mp4', 'MPEG-4');
vid_out_obj.FrameRate = vid_obj.FrameRate;
open(vid_out_obj);
%figure()
i = 1;
while hasFrame(vid_obj)
    frame = readFrame(vid_obj);
    frame = gpuArray(frame);
    frame = rgb2hsv(frame);
    mask = isbetween(frame(:,:,1), hue_bound(1), hue_bound(2)) & ...
             isbetween(frame(:,:,2), sat_bound(1), sat_bound(2)) & ...
             isbetween(frame(:,:,3), val_bound(1), val_bound(2));
    frame(:,:,3) = mask.*frame(:,:,3);
    frame = hsv2rgb(frame);
    frame = gather(frame);
    writeVideo(vid_out_obj,frame);
    %frame(:,:,1) = frame(:,:,1).*mask;
    
    %frame = (frame(:,:,2) + frame(:,:,3))/2;
    %imshow(frame);
    %impixelinfo;

    % Display the progress
    percentDone = 100 * i / vid_obj.NumFrames;
    % Create the message string, e.g., "Percent done: 5.1"
    msg = sprintf('Percent done: %3.1f', percentDone); 
    
    % Print the backspaces first to erase the previous message, 
    % then print the new message
    fprintf([reverseStr, msg]);
    
    % Store the string of backspaces needed to erase the current message next time
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    i = i+1;
end
close(vid_out_obj);
%}