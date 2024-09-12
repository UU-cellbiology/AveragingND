[![](https://github.com/ekatrukha/AveragingND/actions/workflows/build-main.yml/badge.svg)](https://github.com/ekatrukha/AveragingND/actions/workflows/build-main.yml)

AveragingND
===

Averaging of ND (2/3/4D) images using normalized cross-correlation.   
Work in progress, use at your own risk.  

This plugin uses masked normalized cross-correlation calculation implemented according to the paper:  
_D. Padfield, “[Masked object registration in the Fourier domain](https://doi.org/10.1109/TIP.2011.2181402)” IEEE Transactions on Image Processing (2012)._  
So basically it can account (exclude) image values that are zero and therefore can account for arbitrary image dimensions during fast CC calculation using FFT.  

To install it, you need to copy a jar file from the [latest release](https://github.com/ekatrukha/AveragingND/releases) to _jars_ folder of your [FIJI](https://fiji.sc/) installation. The plugin commands should appear in the _Plugins->AveragingND_ menu.   

The plugin works with 2D/3D/4D multichannel images.  

## Iterative averaging

Provided with a set of opened in FIJI or stored on hard drive (.tif) images of the same dimensions, the plugin will perform iterative registration/averaging.  
The images could be of differen sizes, but should have the same dimensions.  
At the initial step it will create an reference image as an average of input images (depending on initial template settings).  
After that it will iteratively register input images to this average, update it, register again.   
Depending on the number of iterations or convergence of average CC it will stop and show the final image.  
The precision of registration is up to one voxel.  
 
## Single template registration

This command registers template image(dataset) to the reference image using rigid translation (in all axis).  
The precision of registration is up to one voxel (subpixel is not available).   
It is possible to limit the maximum shift of the template by a fraction of its size (in all dimensions).  
In case of multi-channel images, only one user selected channel from each image will be used for the registration.  


## Pairwise CC

Given a set of input images, the plugin will pairwise calculate the maximum of cross-correlation among them.  

----------

Developed in <a href='http://cellbiology.science.uu.nl/'>Cell Biology group</a> of Utrecht University.  
<a href="mailto:katpyxa@gmail.com">E-mail</a> for any questions.




