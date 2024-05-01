package tests;


import averagingND.FloatTiffImgWrap;
import averagingND.GenNormCC;
import averagingND.MaskedNormCC;
import averagingND.MiscUtils;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.FinalInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.StopWatch;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class optimizeCC 
{
	public static void main( final String[] args )
	{

		new ImageJ();

		//ImagePlus image1 = IJ.openImage("/home/eugene/Desktop/projects/UnequalTiffs/BB/001ch2.tif");
		//ImagePlus image2 = IJ.openImage("/home/eugene/Desktop/projects/UnequalTiffs/BB/002ch2.tif");
		ImagePlus image1 = IJ.openImage("/home/eugene/Desktop/projects/UnequalTiffs/BB/001ch2.tif");
		ImagePlus image2 = IJ.openImage("/home/eugene/Desktop/projects/UnequalTiffs/BB/002ch2.tif");

		Img<FloatType> ref = ImagePlusAdapter.convertFloat(image1);
		Img<FloatType> temp = ImagePlusAdapter.convertFloat(image2);
		
		StopWatch stopwatch;
		ImageJFunctions.show(ref,"ref");
		ImageJFunctions.show(temp,"temp");
		

		FinalInterval limit = new FinalInterval(new long[] {-20,-20,-20},new long[] {20,20,20});
		
		MaskedNormCC newCC = new MaskedNormCC();
		newCC.bCenteredLimit = true;
		newCC.limInterval = limit;
		newCC.bZeroMask = true;
		newCC.bVerbose = true;		
		stopwatch = StopWatch.createAndStart();		
		newCC.caclulateMaskedNormCC(ref, temp, true);
		System.out.println( "new CC time: " + stopwatch );
		
		GenNormCC oldCC = new GenNormCC();
		oldCC.bCenteredLimit = true;
		oldCC.limInterval = limit;
		oldCC.bZeroMask = true;
		oldCC.bVerbose = true;
		stopwatch= StopWatch.createAndStart();
		oldCC.caclulateGenNormCC(ref, temp, true);
		System.out.println( "old CC time: " + stopwatch );
		
	}

}
