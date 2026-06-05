file = "1500f1";

mov_name = file + ".mov";
% create video reader object for reading teh video files
v = VideoReader(mov_name); 
% Display first frame and ask the user to draw a box covering the lazer dot
frame = readFrame(v); % Read first frame

i = 1;

while hasFrame(v)
    img = readFrame(v);
    filename = sprintf("%03d",i)+".jpg";
    fullname = file + "_frames/" + filename;
    imwrite(img,fullname)    % Write to a JPEG file (001.jpg, 002.jpg, ..., 121.jpg)
    i = i+1;
end