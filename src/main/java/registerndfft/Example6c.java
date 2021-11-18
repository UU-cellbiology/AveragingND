package registerndfft;

import ij.ImageJ;

import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;

import net.imglib2.algorithm.fft.FourierConvolution;
import net.imglib2.algorithm.fft.FourierTransform;
import net.imglib2.algorithm.fft.InverseFourierTransform;
import net.imglib2.converter.ComplexImaginaryFloatConverter;
import net.imglib2.converter.ComplexPhaseFloatConverter;
import net.imglib2.converter.ComplexRealFloatConverter;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.RealSum;

/**
 * Perform template matching by convolution in the Fourier domain
 *
 * @author Stephan Preibisch
 * @author Stephan Saalfeld
 */
public class Example6c
{
	public Example6c() throws ImgIOException, IncompatibleTypeException
	{
		// open with SCIFIO ImgOpener as FloatTypes
		ImgOpener io = new ImgOpener();
		/*final Img< FloatType > image = io.openImgs( "DrosophilaWing.tif",
			new FloatType() ).get( 0 );
		final Img< FloatType > template = io.openImgs( "WingTemplate.tif",
			new FloatType() ).get( 0 );*/
		final Img< FloatType > image = io.openImgs( "oval.tif",
		new FloatType() ).get( 0 );
	final Img< FloatType > template = io.openImgs( "oval_13_15.tif",
		new FloatType() ).get( 0 );

		// display image and template
		ImageJFunctions.show( image ).setTitle( "input" );
		ImageJFunctions.show( template ).setTitle( "template" );

		// compute fourier transform of the template
		final FourierTransform< FloatType, ComplexFloatType > fft =
			new FourierTransform<>(template, new ComplexFloatType() );
		fft.process();
		final Img< ComplexFloatType > templateFFT = fft.getResult();

		// display fft (by default in generalized log power spectrum
		ImageJFunctions.show( templateFFT ).setTitle( "fft power spectrum" );
		// display fft phase spectrum
		ImageJFunctions.show( templateFFT,
			new ComplexPhaseFloatConverter< ComplexFloatType >() )
				.setTitle( "fft phase spectrum" );
		// display fft real values
		ImageJFunctions.show( templateFFT,
			new ComplexRealFloatConverter< ComplexFloatType >() )
				.setTitle( "fft real values" );
		// display fft imaginary values
		ImageJFunctions.show( templateFFT,
			new ComplexImaginaryFloatConverter< ComplexFloatType >() )
				.setTitle( "fft imaginary values" );

		// complex invert the kernel
		final ComplexFloatType c = new ComplexFloatType();
		for ( final ComplexFloatType t : templateFFT )
		{
			c.set( t );
			t.complexConjugate();
			c.mul( t );
			t.div( c );
		}

		// compute inverse fourier transform of the template
		final InverseFourierTransform< FloatType, ComplexFloatType > ifft =
			new InverseFourierTransform<>( templateFFT, fft );
		ifft.process();
		final Img< FloatType > templateInverse = ifft.getResult();

		// display the inverse template
		ImageJFunctions.show( templateInverse ).setTitle( "inverse template" );

		// normalize the inverse template
		//norm( templateInverse );

		// compute fourier convolution of the inverse template and the image and display it
		ImageJFunctions.show( FourierConvolution.convolve( image, templateInverse ) );
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
	 * Norms all image values so that their sum is 1
	 *
	 * @param iterable - the image data
	 */
	public static void norm( final Iterable< FloatType > iterable )
	{
		final double sum = sumImage( iterable );

		for ( final FloatType type : iterable )
			type.setReal( type.get() / sum );
	}

	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();

		// run the example
		new Example6c();
	}
}
