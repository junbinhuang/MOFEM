%% Two end colors:
minrgb=[0.5 0 1];
maxrgb=[1 0 0.5];

%% Convert to HSV colors.
minhsv=rgb_to_hsv(minrgb);
maxhsv=rgb_to_hsv(maxrgb);

%% Number of colors:
numc=102; % Use 102 to match the adina color!!!
t=linspace(0,1,numc);

global myColor;

%% My colormap:
myColor=zeros(numc-2,3);

for i=1:numc-2
    myColor(i,:)=hsv_to_rgb(t(i+1)*maxhsv+(1-t(i+1))*minhsv);
end

clear t i numc minrgb maxrgb minhsv maxhsv

% colormap(myColor);

%% Two functions used here:
function [hsv] = rgb_to_hsv(rgb) %see Foley & Van Dam Fundamentals of interactive computer graphics P. 615
    r=rgb(1);
    g=rgb(2);
    b=rgb(3);
    maxc=max(rgb);
    minc=min(rgb);
    v=maxc;
    if maxc~=0
        s=(maxc-minc)/maxc;
    else
        s=0;
    end
    if s==0
        h=0;
    else
        rc=(maxc-r)/(maxc-minc);
        gc=(maxc-g)/(maxc-minc);
        bc=(maxc-b)/(maxc-minc);
        if r==maxc
            h=bc-gc;
        else
            if g==maxc
                h=2+rc-bc;
            else
                if b==maxc
                    h=4+gc-rc;
                end
            end
        end
    end
    h=h*60;
%     if h<0 % Moved to the next function.
%         h=h+360;
%     end
    hsv=[h s v];
end

function [rgb] = hsv_to_rgb(hsv) %see Foley & Van Dam Fundamentals of interactive computer graphics P. 616
    h=hsv(1);
    s=hsv(2);
    v=hsv(3);
    if h<0
        h=h+360;
    end
    if s==0
        r=v;
        g=v;
        b=v;
    else
        if h==360
            h=0;
        end
        h=h/60;
        i=floor(h);
        f=h-i;
        p=v*(1-s);
        q=v*(1-s*f);
        t=v*(1-s*(1-f));
        switch i
            case 0
                r=v;
                g=t;
                b=p;
            case 1
                r=q;
                g=v;
                b=p;
            case 2
                r=p;
                g=v;
                b=t;
            case 3
                r=p;
                g=q;
                b=v;
            case 4
                r=t;
                g=p;
                b=v;
            case 5
                r=v;
                g=p;
                b=q;
        end
    end
    rgb=[r g b];
end