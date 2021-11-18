package registerndfft;

import ij.IJ;
import ij.ImageJ;

import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;

import net.imglib2.algorithm.integral.IntegralImgDouble;
import net.imglib2.algorithm.fft2.FFT;
import net.imglib2.algorithm.fft2.FFTMethods;
import net.imglib2.Cursor;
//import net.imglib2.Dimensions;

import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
//import net.imglib2.algorithm.fft.FourierConvolution;
//import net.imglib2.algorithm.fft.FourierTransform;
//import net.imglib2.algorithm.fft.InverseFourierTransform;
import net.imglib2.converter.ComplexImaginaryFloatConverter;
import net.imglib2.converter.ComplexPhaseFloatConverter;
import net.imglib2.converter.ComplexRealFloatConverter;
import net.imglib2.converter.RealDoubleConverter;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.RealSum;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import net.imglib2.view.composite.RealComposite;

/**
 * Perform template matching by convolution in the Fourier domain
 *
 * @author Stephan Preibisch
 * @author Stephan Saalfeld
 */
public class Example_fft2
{
	public Example_fft2() throws ImgIOException, IncompatibleTypeException
	{
		
		int i;
		
		// open with SCIFIO ImgOpener as FloatTypes
		ImgOpener io = new ImgOpener();
		
		/*final Img< FloatType > image = io.openImgs( "Zstack 1.1-1-1.tif",
			new FloatType() ).get( 0 );
		final Img< FloatType > template = io.openImgs( "Zstack 1.2-2-1.tif",
			new FloatType() ).get( 0 );
			*/
			
		
		final Img< FloatType > image = io.openImgs( "linetest2D.tif",
				new FloatType() ).get( 0 );
		final Img< FloatType > template = io.openImgs( "linetest2D_shift.tif",
				new FloatType() ).get( 0 );

		int nDim = image.numDimensions();
		
		long [] imgDim = new long[nDim];
		long [] temDim = new long[nDim];

		image.dimensions(imgDim);
		template.dimensions(temDim);
		
		long nPixImg = imgDim[0];
		long nPixTem = temDim[0];
		
		for(i=1;i<nDim;i++)
		{
			nPixImg *= imgDim[i];
			nPixTem *= temDim[i];			
		}
		
		
		//subtract average values
		subtract_mean(image,nPixImg);
		subtract_mean(template,nPixTem);
		
		// display image and template
		ImageJFunctions.show( image ).setTitle( "input" );
		ImageJFunctions.show( template ).setTitle( "template" );
		
		if(image.numDimensions()!=template.numDimensions())
		{
			IJ.log("different dimensions of input and template!");
			return;
		}

		long [][] cropCorr = new long[2][nDim];
		long [][] finDim = new long[2][nDim];
		for(i=0;i<nDim;i++)
		{
			cropCorr[0][i]=Math.round((-0.5)*temDim[i]);
			cropCorr[1][i]=Math.round(0.5*imgDim[i]);
			finDim[0][i]=0;
			finDim[1][i]=Math.max(imgDim[i], temDim[i]);
		}
		
		
		
		
		//long[] paddedDimensions = new long[template.numDimensions()];
		//long[] fftSize = new long[template.numDimensions()];
		
		//FFTMethods.dimensionsRealToComplexFast(template, paddedDimensions, fftSize);
		
		final ImgFactory< ComplexFloatType > factoryComplex = new ArrayImgFactory< ComplexFloatType >();
		final ImgFactory< FloatType > factoryFloat = new ArrayImgFactory< FloatType >();
		
		final Img< ComplexFloatType > imageFFT2;
		final Img< ComplexFloatType > templateFFT2;
		FinalInterval padDim = new FinalInterval( finDim[0] ,  finDim[1] );
		//padd both images to zero
		IntervalView< FloatType > padIm = Views.interval(Views.extendZero(image),padDim);
		IntervalView< FloatType > padTem = Views.interval(Views.extendZero(template),padDim);
		
		imageFFT2=FFT.realToComplex(padIm, factoryComplex);
		templateFFT2=FFT.realToComplex(padTem, factoryComplex);
		
		//calculate sum of original FFT and conjugate FFT of the template 
		FFTMethods.complexConjugate(templateFFT2);		
		final Cursor< ComplexFloatType > imC = imageFFT2.cursor();
		final Cursor< ComplexFloatType > templateC = templateFFT2.cursor();		
		final ComplexFloatType cTemp = new ComplexFloatType();
		while(imC.hasNext())
		{
			imC.fwd();
			templateC.fwd();
			cTemp.set(imC.get());
			cTemp.mul(templateC.get());
			templateC.get().set(cTemp);
		}
				
		//return back to normal space
		final Img< FloatType > invertedT = FFT.complexToReal(templateFFT2, factoryFloat, new FloatType());		
		
		//range of CC
		FinalInterval interval = new FinalInterval( cropCorr[0] ,  cropCorr[1] );

		//swap quadrants
		IntervalView< FloatType > ivCC = Views.interval( Views.extendPeriodic( invertedT ),interval);
		
		//calculate summed-area tables for padded images (intensity squared)
		square_values(padIm);
		square_values(padTem);
		IntegralImgDouble<FloatType> integralImgProc = new IntegralImgDouble<FloatType>(padIm,new DoubleType(),new RealDoubleConverter<FloatType>());
		IntegralImgDouble<FloatType> integralTemProc = new IntegralImgDouble<FloatType>(padTem,new DoubleType(),new RealDoubleConverter<FloatType>());
		
		RandomAccessibleInterval<DoubleType> imgIntegral = null;
		RandomAccessibleInterval<DoubleType> temIntegral = null;
		if(integralImgProc.process())
		{
			imgIntegral = integralImgProc.getResult();
		}
		if(integralTemProc.process())
		{
			temIntegral = integralTemProc.getResult();
		}
		
		/*// let's take a look at one of them
		ImageJFunctions.show(imgIntegral).setTitle( "integral_image" );
		*/
		
		
		
		Point shift = new Point(nDim);
		computeMaxLocation(ivCC,shift);
		for(i=0;i<nDim;i++)
		{
			IJ.log("dim "+Integer.toString(i)+": "+Integer.toString(shift.getIntPosition(i)));
		}
		Views.interval(Views.extendZero(template),interval);
		long [] shift_pos=new long[nDim];
		shift.localize(shift_pos);
		Views.interval(Views.translate(Views.extendZero(template), shift_pos),template);
		
		ImageJFunctions.show(ivCC).setTitle( "cross-correlation" );
		ImageJFunctions.show(Views.interval(Views.translate(Views.extendZero(template), shift_pos),template)).setTitle("registered template");

	}
	
	public < T extends Comparable< T > & Type< T > > void computeMaxLocation(
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
	 * Norms all image values so that their sum is 1
	 *
	 * @param iterable - the image data
	 */
	public static void subtract_mean( final Iterable< FloatType > iterable, long nPixelsN)
	{
		final double sum = sumImage( iterable );
		final double meanV =sum/(double)nPixelsN; 
		for ( final FloatType type : iterable )
			type.setReal( type.get() - meanV );
	}
	/**
	 * Norms all image values so that their sum is 1
	 *
	 * @param iterable - the image data
	 */
	public static void square_values( final Iterable< FloatType > iterable)
	{
		float val=0;
		for ( final FloatType type : iterable )
		{
			val = type.get();
			type.setReal( val*val );
		}
	}
		
	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();

		// run the example
		new Example_fft2();
	}
}
