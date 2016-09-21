function [l,c] = image_coordinates(h,i)

l = mod(i,h);

if l==0 
    l=h; 
end

c = (i - l)/h + 1;
    