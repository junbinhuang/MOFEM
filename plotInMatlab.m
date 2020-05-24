%% Read data from the Python outputs.
fid=fopen('nSampling.txt');
data = textscan(fid, '%d', 'CommentStyle','#', 'CollectOutput',true);
fclose(fid);
numbers=cell2mat(data);

fid=fopen('boundaryNumber.txt');
data = textscan(fid, '%d', 'CommentStyle','#', 'CollectOutput',true);
fclose(fid);
boundaryNumber=cell2mat(data);

fid=fopen('boundaryCoord.txt');
data = textscan(fid, '%f %f', 'CommentStyle','#', 'CollectOutput',true);
fclose(fid);
bCoord=cell2mat(data);

nSampling=numbers(1); nComponent=numbers(2);

a = '%f ';
for i = 1:nSampling-2 % Assume nSampling is at least 2.
    a = [a,'%f '];
end
a = [a,'%f'];
fid=fopen('SampleData.txt');
data = textscan(fid, a, 'CommentStyle','#', 'CollectOutput',true);
fclose(fid);
totalData=cell2mat(data);

clear data a

%% Now we can start plotting figures.
run ADINAColor
% Plot the boundary.
for iFigure=1:nComponent
    %% Plot only the effective stress.
%     if iFigure~=6
%         continue
%     end
    %%
    
    figure
    hold on
    axis equal
    axis off

    boundaryCoord=bCoord;

    %% Plot the boundary.
    for i=1:length(boundaryNumber)
        coord=boundaryCoord(1:boundaryNumber(i),:);
        boundaryCoord=boundaryCoord(boundaryNumber(i)+1:end,:);

        plot(coord(:,1),coord(:,2),'k','LineWidth',3)
        if coord(1,1)~=coord(end,1) || coord(1,2)~=coord(end,2)
            plot(coord([1,end],1),coord([1,end],2),'k','LineWidth',3)
        end
    end

    % Plot the contourf results!
    index=getIndex(iFigure,nSampling,nComponent,size(totalData,1));
    plotData=totalData(index,:);
    
    nScale=100;
    nDomain=size(plotData,1)/nSampling/3;

    xArray=plotData(1:3:end,:);
    yArray=plotData(2:3:end,:);
    zArray=plotData(3:3:end,:);

    xMin=min(xArray(:));xMax=max(xArray(:));
    yMin=min(yArray(:));yMax=max(yArray(:));
    zMin=min(zArray(:));zMax=max(zArray(:));
    
    %% Set limits
%     if iFigure==5
%         zMin=-20000;
%         zMax=20000;
%     end
%     if iFigure==2
%         zMin=-3e-4;
%         zMax=0;
%     end
%     if iFigure==6
%         zMin=0;
%         zMax=5000;
%     end
    %%

    scale=linspace(zMin,zMax,nScale);

    for i=1:nDomain
%         if min(xArray(nSampling*(i-1)+1:nSampling*i,:))<23 % For the
%         % stress inside inclusion.
%             continue
%         end
%         if max(xArray(nSampling*(i-1)+1:nSampling*i,:))>31
%             continue
%         end
%         if min(yArray(nSampling*(i-1)+1:nSampling*i,:))<-4
%             continue
%         end
%         if max(yArray(nSampling*(i-1)+1:nSampling*i,:))>4
%             continue
%         end
        
        myContourf(xArray(nSampling*(i-1)+1:nSampling*i,:),...
                   yArray(nSampling*(i-1)+1:nSampling*i,:),...
                   zArray(nSampling*(i-1)+1:nSampling*i,:),scale)
    end

    xlim([xMin,xMax])
    ylim([yMin,yMax])
%     xlim([23.5,30.5])
%     ylim([-3.5,3.5])
    %% plot the circle.
%     theta=linspace(0,2*pi,100);
%     coord=zeros(100,2);
%     for i=1:100
%         coord(i,1)=27+3*cos(theta(i));
%         coord(i,2)=3*sin(theta(i));
%     end
%     plot(coord(:,1),coord(:,2),'k','LineWidth',0.8)
end

clear plotData xArray yArray zArray scale nDomain xMin xMax yMin yMax zMin zMax...
    iFigure i fid coord

%% Some functions used:
function myContourf(x,y,z,scale)
%Used in visualization
    global myColor;
    contourf(x,y,z,scale,'LineStyle','none');
    set(gca,'ticklabelinterpreter','latex','fontsize',21)
    colormap(myColor);
    colorbar('ticklabelinterpreter','latex')
end

function index=getIndex(i,nSampling,nComponent,n)
    index=zeros(1,n/nComponent);
    
    for j=1:n/nComponent/3/nSampling
        index(1+(j-1)*3*nSampling:3*j*nSampling)=(1:3*nSampling)+...
            (i-1)*3*nSampling+(j-1)*nComponent*3*nSampling;
    end
end