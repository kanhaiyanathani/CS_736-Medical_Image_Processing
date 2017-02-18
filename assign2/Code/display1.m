function []= display1(im)
% im=(im-min(im(:)))/(max(im(:))-min(im(:)));
im(im<0)=0;
im = im/max(im(:));
imshow(im);
end