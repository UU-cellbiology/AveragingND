package averagingND;


import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class testPeriodic {
	/**
	 * Show 16-bit volume.
	 */
	public static void main( final String[] args )
	{
		new ImageJ();
		final ImagePlus imp = IJ.openImage( "https://imagej.nih.gov/ij/images/blobs.gif" );
		final Img< FloatType > img = ImageJFunctions.convertFloat( imp );
		ImageJFunctions.show(img).setTitle("Original");
		double degree = 0.7;
		
		int nDim = img.numDimensions();
		long [] dims = new long [nDim];
		img.dimensions(dims);
		
		long [][] newBounds = new long [2][nDim];
		for(int i=0;i<nDim;i++)
		{
			newBounds[0][i]= Math.round((-1)*degree*dims[i]);
			newBounds[1][i]= Math.round(degree*dims[i]);
		}
		FinalInterval interval = new FinalInterval( newBounds[0] ,  newBounds[1] );
		
		//let's do periodic 
		IntervalView< FloatType > ivPeriod = Views.interval( Views.extendPeriodic( img ),interval);
		//show it
		ImageJFunctions.show(ivPeriod).setTitle("Periodic");
		
		Cursor <FloatType> ivCursor = ivPeriod.cursor();
		float newVal=0.0f;
		while (ivCursor.hasNext())
		{
			ivCursor.fwd();
			//invert
			newVal = ivCursor.get().get();
			newVal=255-newVal;
			ivCursor.get().set(newVal);
				
		}
		ImageJFunctions.show(ivPeriod).setTitle("Periodic_inverted");
	}
}
