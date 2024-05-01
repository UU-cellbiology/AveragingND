package tests;

import averagingND.FloatTiffImgWrap;
import averagingND.MiscUtils;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;


public class virtualWrapTest 
{
	public static void main( final String[] args )
	{
		
		new ImageJ();
		ImagePlus imageIn = IJ.openImage("/home/eugene/Desktop/projects/UnequalTiffs/BB/001f.tif");
		Img< FloatType >img = FloatTiffImgWrap.wrapVirtualFloat(imageIn, MiscUtils.getDimensionsTextImageJ(imageIn));
		ImageJFunctions.show(img, "test_raw");
		
	}
}
