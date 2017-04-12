
# Image augmentation in R

In deep learning, we often employ data augmentation techniques to increase sample size. Data augmentation for image include flip, flop, random rotation, crop and ZCA whitening etc. In this brief tutorial, I am using R package [`OpenImageR`](https://cran.r-project.org/package=OpenImageR) to perform image augmentation. The following example assumes that you have successfully installed `OpenImageR`.

## Install `OpenImageR`
```R
>install.package("OpenImageR")
>library(OpenImageR) # loading install package
```

## An R routine to perform image augmentation
```R
image_augmentation <- function(filename) {

    outfn.prefix=sub(".jpg$|.jpeg$|.png$|.tiff$","",fl, ignore.case = TRUE)

    image=readImage(fl)

    horizontal_flipped_image=flipImage(image, mode = "horizontal")
    writeImage(horizontal_flipped_image, paste(outfn.prefix, "_horizontal_flipped.jpg",sep=""))

    vertical_flipped_image=flipImage(image, mode = "vertical")
    writeImage(vertical_flipped_image, paste(outfn.prefix, "_vertical_flipped.jpg",sep=""))

    # randomly rotation
    rotated_image=rotateImage(image, angle = sample(1:359,1))
    writeImage(rotated_image, paste(outfn.prefix, "_rotated.jpg",sep=""))

    # image cropping, the cropping dimension can be changed.
    cropped_image=cropImage(image, dim(image)[1]*0.62, dim(image)[2]*0.62)
    writeImage(cropped_image, paste(outfn.prefix, "_cropped.jpg",sep=""))
	
    # ZCA whitening
    zcawhitening_image=ZCAwhiten(image, k = 100, epsilon = 0.1)
    writeImage(zcawhitening_image, paste(outfn.prefix, "_zcawhitening.jpg",sep=""))

    print(filename)
} 
```

## Perform image augmentation in parallel
```R
library(parallel)
library(OpenImageR)

# `fn` contains path to every image file, one per line
fls=read.table(fn,header=FALSE,stringsAsFactors=FALSE)$V1

# r=mclapply(fls, function(fl) {
	image_augmentation(fl)
}, mc.cores=60)
```


# Processing DICOM formatted medical images
In this brief tutorial, the Matlab Image Processing Toolbox is used to process medical image in DICOM format.

```matlab
[X, map] = dicomread('~/tmp.dcm'); % read dicom file
imshow(X)
```

