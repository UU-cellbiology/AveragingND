package fiji.plugin.RegisterNDFFT;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.img.Img;
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
	public static < T extends Comparable< T > & Type< T > > void computeMaxLocation(
			final IterableInterval< T > input, final Point maxLocation )
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
}
