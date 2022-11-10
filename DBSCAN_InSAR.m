function [class, lonlat_index] = DBSCAN_InSAR(velocity, lonlat, MinPts,Eps,...
    velocity_thre,lon_lim,lat_lim)

%% Density-based spatial clustering of applications with noise (DBSCAN) for InSAR, edited by Chuang Song
%% references: Ester M, Kriegel H-P, Sander J, Xu X. A density-based algorithm for discovering clusters in large spatial databases with noise. In: Proc. 2nd Int. Conf. Knowledge Discovery and Data Mining (KDD’96)). AAAI Press (1996).
%% references: https://cloud.tencent.com/developer/article/1064264, credited by zhaozhiyong

% Inputs -->
% velocity: vectors of InSAR velocity
% lonlat:  vector of [longitude latitude]
% MinPts: minimum number of pixels in the searching window, e.g., 3
% Eps: radius (m) of the searching window, e.g., 100

% Optional --> 
% velocity_thre: threshold of InSAR velocity (mm/yr) for determining pixel similarity, default 9999 to neglect this constraint.
% lon_lim: [minlon maxlon], lat_lim: [minlat, maxlat]: constraint of spatial range

% Outputs -->
% class: with nans
% lonlat_index: just the index for spatial range and notnan velocity

% after this script -->
% >> to extract valid classes, firstly, class_index = find(~isnan(class)); 
% >> when without lon_lim and lat_lim: lonlat_output = lonlat(class_index,:)
% >> when with lon_lim and lat_lim: lonlat_inside = lonlat(lonlat_index,:); lonlat_output = lonlat_inside(class_index,:);

%% Import data
v_mm_yr = velocity;
lon0 = lonlat(1,1); lat0 = lonlat(1,2);
xy = llh2local(lonlat',[lon0,lat0])'*1000; % m

if nargin > 4
    velocity_thre = 9999;
end
lon_in = [];
lat_in = [];
if nargin > 5
    lon_in = lon_lim;
end
if nargin > 6
    lat_in = lat_lim;
end
if ~isempty(lon_in)
    lon_index = lonlat(:,1)>=lon_in(1)&lonlat(:,1)<=lon_in(2);
else
    lon_index = 1:length(v_mm_yr);
end
if ~isempty(lat_in)
    lat_index = lonlat(:,2)>=lat_in(1)&lonlat(:,2)<=lat_in(2);
else
    lat_index = 1:length(v_mm_yr);
end
notnan_index = ~isnan(velocity);
lonlat_index = lon_index & lat_index & notnan_index';


data = xy(lonlat_index,:);
dataz = v_mm_yr(lonlat_index);

[m,~] = size(data);

x = [(1:m)' data];
[m,n] = size(x);
types = zeros(1,m);% 1: Core point 1, 0: boundary point, -1: noise point
dealed = zeros(m,1);% 0: not processed, 1: processed
class = -9999*ones(m,1);

fprintf('Calculate distance matrix for %d pixels ...\n',length(dataz));
dis = calDistance(x(:,2:n));
distance_matrix = dis;
number = 1;% to mark different classes

fprintf('Begin clustering ...\n');
%% pixel by pixel
for i = 1:m
    if mod(i,1000) == 0
        fprintf('Clustering No.%d\n',i);
    end
    % find un-processed pixels
    if dealed(i) == 0
        xTemp = x(i,:);
        D = dis(i,:);
        ind0 = find(D<=Eps);% find all pixels inside Eps
        ind = ind0(abs(dataz(i)-dataz(ind0))<=velocity_thre);
        
        %% identify the type of pixels
        
        % Boundary point
        if length(ind) > 1 && length(ind) < MinPts
            types(i) = 0;
            class(i) = 0;
        end
        % Noise point
        if length(ind) == 1
            types(i) = -1;
            class(i) = -1;
            dealed(i) = 1;
        end
        % Core point
        if length(ind) >= MinPts
            types(xTemp(1,1)) = 1;
            class(ind) = number;
            
            % Judge if density-reachable
            while ~isempty(ind)
                yTemp = x(ind(1),:);
                dealed(ind(1)) = 1;
                ind(1) = [];
                
                ind = unique(ind);
                
                D = dis(yTemp(1,1),:);
                ind_1_0 = find(D<=Eps);
                ind_1 = ind_1_0(abs(dataz(yTemp(1,1))-dataz(ind_1_0))<=velocity_thre);
                
                if length(ind_1)>1 % process non-noise points
                    class(ind_1) = number;
                    if length(ind_1) >= MinPts
                        types(yTemp(1,1)) = 1;
                    else
                        types(yTemp(1,1)) = 0;
                    end
                    
                    for j=1:length(ind_1)
                        if dealed(ind_1(j)) == 0
                            dealed(ind_1(j)) = 1;
                            ind=[ind ind_1(j)];
                            class(ind_1(j))=number;
                        end
                    end
                end
            end
            number = number + 1;
        end
    end
end

% set remaining un-processed points as noise points
ind_2 = find(class==0);
class(ind_2) = -1;
types(ind_2) = -1;

%% Visualisation
% % % % plotscatter(lonlat(lonlat_index,1),lonlat(lonlat_index,2),dataz,[-100 100],0,0,8);
class(class==-1) = nan;
% % % % plotscatter(lonlat(lonlat_index,1),lonlat(lonlat_index,2),class,[0 max(class)],0,0,8);
