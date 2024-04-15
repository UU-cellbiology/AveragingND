package averagingND;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import ij.ImagePlus;
import ij.measure.Calibration;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.RealSum;

public class MiscUtils {
	
	/**
	 * Generic, type-agnostic method to create an identical copy of an Img
	 *
	 * @param input - the Img to copy
	 * @return - the copy of the Img
	 */
	public static < T extends Type< T > > Img< T > copyImage( final Img< T > input )
	{
		// create a new Image with the same properties
		// note that the input provides the size for the new image as it implements
		// the Interval interface
		Img< T > output = input.factory().create( input );
 
		// create a cursor for both images
		Cursor< T > cursorInput = input.cursor();
		Cursor< T > cursorOutput = output.cursor();
 
		// iterate over the input
		while ( cursorInput.hasNext())
		{
			// move both cursors forward by one pixel
			cursorInput.fwd();
			cursorOutput.fwd();
 
			// set the value of this pixel of the output image to the same as the input,
			// every Type supports T.set( T type )
			cursorOutput.get().set( cursorInput.get() );
		}
 
		// return the copy
		return output;
	}
	
	
	/** computes location of the maximum voxel **/
	public static < T extends Comparable< T > & Type< T > > T computeMaxLocation(
			final IterableInterval< T > input, final Point maxLocation)
		{
			// create a cursor for the image (the order does not matter)
			final Cursor< T > cursor = input.cursor();
	 
			// initialize min and max with the first image value
			T type = cursor.next();
			T max = type.copy();
	 
			// loop over the rest of the data and determine min and max value
			while ( cursor.hasNext() )
			{
				// we need this type more than once
				type = cursor.next();
	 
	 
				if ( type.compareTo( max ) > 0 )
				{
					max.set( type );
					maxLocation.setPosition( cursor );
				}
			}
			return max;
		}
	
	/**
	 * Computes the sum of all pixels in an iterable using RealSum
	 *
	 * @param iterable - the image data
	 * @return - the sum of values
	 */
	public static < T extends RealType< T > > double sumImage( final Iterable< T > iterable )
	{
		final RealSum sum = new RealSum();

		for ( final T type : iterable )
			sum.add( type.getRealDouble() );

		return sum.getSum();
	}
	
	/**
	 * subtracts mean valut from the image 
	 *
	 * @param iterable - the image data
	 * @param nPixelsN - number of pixels in the image
	 */
	public static void subtractMean( final Iterable< FloatType > iterable, long nPixelsN)
	{
		final double sum = sumImage( iterable );
		final double meanV =sum/(double)nPixelsN; 
		for ( final FloatType type : iterable )
			type.setReal( type.get() - meanV );
	}
	
	/**
	 * squares all image values 	 
	 * @param iterable - the image data
	 */
	public static void squareValues( final Iterable< FloatType > iterable)
	{
		float val=0;
		for ( final FloatType type : iterable )
		{
			val = type.get();
			type.setReal( val*val );
		}
	}
	
	/** assume there is no time axis **/
	public static ImagePlus wrapFloatImgCal(RandomAccessibleInterval< FloatType > img, String sTitle, Calibration cal, boolean nCh)
	{
		ImagePlus outIP = ImageJFunctions.wrap(img,sTitle);
		
		
		if(img.numDimensions()==3)
		{
			if(nCh)
			{
				outIP.setDimensions((int)img.dimension(2), 1 , 1);
			}
			else
			{
				outIP.setDimensions(1, (int)img.dimension(2), 1);
			}
		}
		if(img.numDimensions()==4)
		{
			outIP.setDimensions((int)img.dimension(2), (int)img.dimension(3), 1);
		}
		if(cal!=null)
		{
			outIP.setCalibration(cal);
		}
		return outIP;
		//outIP.show();
	}
	
    public static List<String> findFiles(Path path, String fileExtension)
            throws IOException {

            if (!Files.isDirectory(path)) {
                throw new IllegalArgumentException("Path must be a directory!");
            }

            List<String> result;

            try (Stream<Path> stream = Files.list(path)) {
                result = stream
                        .filter(p -> !Files.isDirectory(p))
                        // this is a path, not string,
                        // this only test if path end with a certain path
                        //.filter(p -> p.endsWith(fileExtension))
                        // convert path to string first
                        .map(p -> p.toString())
                        .filter(f -> f.endsWith(fileExtension))
                        .collect(Collectors.toList());
            }

            return result;
    }
	
    /** returns dimensions of the ImagePlus in XYZTC format**/
	public static String getDimensionsText(final ImagePlus ip)
	{
		String sDims = "XY";
		
		if(ip.getNSlices()>1)
		{
			sDims = sDims + "Z";
		}
		
		if(ip.getNFrames()>1)
		{
			sDims = sDims + "T";
		}
		
		if(ip.getNChannels()>1)
		{
			sDims = sDims + "C";
		}

		return sDims;
	}
}
