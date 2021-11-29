package fiji.plugin.RegisterNDFFT;

import ij.IJ;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.FinalDimensions;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.algorithm.fft2.FFT;
import net.imglib2.algorithm.fft2.FFTMethods;
import net.imglib2.converter.Converters;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.FloatArray;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class GenNormCC {
	
	
	/** dimensionality of images **/	
	//int nDim;
	
	/**
	 * @param image
	 * @param template
	 */
	public static void caclulateGenNormCC(final Img< FloatType > image, final Img< FloatType > template, double max_fraction_shift) //throws ImgIOException, IncompatibleTypeException
	{
		int i;
		int nDim;
		
		
		if(image.numDimensions()!=template.numDimensions())
		{
			IJ.log("different dimensions of input and template!");
			return;
		}
		
		nDim = image.numDimensions();
		long [] imgDim = new long[nDim];
		long [] temDim = new long[nDim];
		image.dimensions(imgDim);
		template.dimensions(temDim);
		final long [] finDim = new long [nDim];
		for(i=0;i<nDim;i++)
		{
			finDim[i]=imgDim[i]+temDim[i]-1;
		}
		
		
		long[] paddedDimensions = new long[template.numDimensions()];
		long[] fftSize = new long[template.numDimensions()];
		FFTMethods.dimensionsRealToComplexFast(new FinalDimensions(finDim), paddedDimensions, fftSize);
		
		System.out.println("done");
		
		//padded intervals
		FinalInterval imgIntPad = (FinalInterval) FFTMethods.paddingIntervalCentered(image, new FinalDimensions(paddedDimensions));
		FinalInterval temIntPad = (FinalInterval) FFTMethods.paddingIntervalCentered(template, new FinalDimensions(paddedDimensions));
		
		//padded versions
		IntervalView< FloatType > padImg = Views.interval(Views.extendZero(image),imgIntPad);
		IntervalView< FloatType > padTem = Views.interval(Views.extendZero(template),temIntPad);
		//ImageJFunctions.show(padImg).setTitle( "padded1" );
		//ImageJFunctions.show(padTem).setTitle( "padded2" );
		
		//unity images
		ArrayImg<FloatType, FloatArray> unityImg = ArrayImgs.floats(imgDim);
		ArrayImg<FloatType, FloatArray> unityTem = ArrayImgs.floats(temDim);
		unityImg.forEach(t->t.set(1.0f));
		unityTem.forEach(t->t.set(1.0f));

		//padded unity images
		IntervalView< FloatType > padUnitImg = Views.interval(Views.extendZero(unityImg),imgIntPad);
		IntervalView< FloatType > padUnitTem = Views.interval(Views.extendZero(unityTem),temIntPad);
		//ImageJFunctions.show(padUnitImg).setTitle( "unipadded1" );
		//ImageJFunctions.show(padUnitTem).setTitle( "unipadded2" );
	
		//squared image values
		Img< FloatType > sqImg = image.copy();
		Img< FloatType > sqTem = template.copy();
		
		sqImg.forEach(t-> t.mul(t.get()));
		sqTem.forEach(t-> t.mul(t.get()));
		//padded squared images
		IntervalView< FloatType > padSqImg = Views.interval(Views.extendZero(sqImg),imgIntPad);
		IntervalView< FloatType > padSqTem = Views.interval(Views.extendZero(sqTem),temIntPad);
		
		//ImageJFunctions.show(padSqImg).setTitle( "sqpadded1" );
		//ImageJFunctions.show(padSqTem).setTitle( "sqpadded2" );
		
		final ImgFactory< ComplexFloatType > factoryComplex = new ArrayImgFactory< ComplexFloatType >(new ComplexFloatType());
		final ImgFactory< FloatType > factoryFloat = new ArrayImgFactory< FloatType >(new FloatType());
		
		//start with FFT
		final Img< ComplexFloatType > imageFFT2 = FFT.realToComplex(padImg, factoryComplex);
		final Img< ComplexFloatType > templateFFT2=FFT.realToComplex(padTem, factoryComplex);
		final Img< ComplexFloatType > unImageFFT2 = FFT.realToComplex(padUnitImg, factoryComplex);
		final Img< ComplexFloatType > unTemplateFFT2=FFT.realToComplex(padUnitTem, factoryComplex);
		final Img< ComplexFloatType > sqImageFFT2 = FFT.realToComplex(padSqImg, factoryComplex);
		final Img< ComplexFloatType > sqTemplateFFT2=FFT.realToComplex(padSqTem, factoryComplex);
		
		//conjugates
		FFTMethods.complexConjugate(templateFFT2);	
		FFTMethods.complexConjugate(unTemplateFFT2);	
		FFTMethods.complexConjugate(sqTemplateFFT2);	
		
		//multiplications
		final Img< ComplexFloatType > deNOM = multCompl(unImageFFT2,unTemplateFFT2);
	    final Img< ComplexFloatType > F1F2 = multCompl(imageFFT2,templateFFT2);
	    final Img< ComplexFloatType > F1I2 = multCompl(imageFFT2,unTemplateFFT2);
	    final Img< ComplexFloatType > I1F2 = multCompl(unImageFFT2,templateFFT2);
	    final Img< ComplexFloatType > F1F1I2 = multCompl(sqImageFFT2,unTemplateFFT2);
	    final Img< ComplexFloatType > I1F2F2 = multCompl(unImageFFT2,sqTemplateFFT2);
	    
	    //inverse FFT
		final Img< FloatType > invdeNOM = FFT.complexToReal(deNOM, factoryFloat, new FloatType());	
		final Img< FloatType > invF1F2 = FFT.complexToReal(F1F2, factoryFloat, new FloatType());	
		final Img< FloatType > invF1I2 = FFT.complexToReal(F1I2, factoryFloat, new FloatType());	
		final Img< FloatType > invI1F2 = FFT.complexToReal(I1F2, factoryFloat, new FloatType());	
		final Img< FloatType > invF1F1I2 = FFT.complexToReal(F1F1I2, factoryFloat, new FloatType());	
		final Img< FloatType > invI1F2F2 = FFT.complexToReal(I1F2F2, factoryFloat, new FloatType());
		
		//ImageJFunctions.show(invdeNOM).setTitle( "denom" );
		//ImageJFunctions.show(invF1F2).setTitle( "invF1F2" );
		//ImageJFunctions.show(invF1I2).setTitle( "invF1I2" );
		//ImageJFunctions.show(invI1F2).setTitle( "invI1F2" );
		//ImageJFunctions.show(invF1F1I2).setTitle( "invF1F1I2" );
		//ImageJFunctions.show(invI1F2F2).setTitle( "invI1F2F2" );
		
		calcTermNum(invF1F2,invF1I2,invI1F2,invdeNOM);
		calcTermDenom(invF1F1I2,invF1I2,invF1I2,invdeNOM);
		calcTermDenom(invI1F2F2,invI1F2,invI1F2,invdeNOM);
		
		//ImageJFunctions.show(invF1F2).setTitle( "term0" );
		//ImageJFunctions.show(invF1F1I2).setTitle( "term1" );
		//ImageJFunctions.show(invI1F2F2).setTitle( "term2" );

		/*
		final Cursor< FloatType > cursFl1 = invF1F2.cursor();
		final Cursor< FloatType > cursFl2 = invF1F1I2.cursor();
		final Cursor< FloatType > cursFl3 = invI1F2F2.cursor();

		double tol = 1E-5;
		while(cursFl1.hasNext())
		{
			cursFl1.fwd();
			cursFl2.fwd();
			cursFl3.fwd();
			
			t2=cursFl2.get().get();
			t3=cursFl3.get().get();
			t2= (float)Math.sqrt(t2*t3);
			if(t2>tol)
			{
				cursFl1.get().mul(1.0f/t2);
			}
			else
			{
				cursFl1.get().set(0.0f);
			}
			//debug
			//cursFl2.get().set(t2);
		}
		*/
		
		
		final Cursor< FloatType > denom1 = invF1F1I2.cursor();
		final Cursor< FloatType > denom2 = invI1F2F2.cursor();
		float maxVal = Float.MIN_VALUE;
		float t2, t3;
		
		
		//calculate denominator and its max value
		//to estimate precision later
		while(denom1.hasNext())
		{
			denom1.fwd();
			denom2.fwd();
			t2 = denom1.get().get();
			t3 = denom2.get().get();
			t2 = (float)Math.sqrt(t2*t3);
			denom1.get().set(t2);
			if(t2>maxVal)
			{
				maxVal=t2;
			}
		}
		
		//numerator and denominator
		double tol = Math.abs(maxVal)*1E-4;//float precision limit
		final Cursor< FloatType > numCurs = invF1F2.cursor();
		final Cursor< FloatType > denomCurs = invF1F1I2.cursor();
		while(numCurs.hasNext())
		{
			numCurs.fwd();
			denomCurs.fwd();
			
			t2 = denomCurs.get().get();
			
			if(t2>tol)
			{
				numCurs.get().mul(1.0f/t2);
			}
			else
			{
				numCurs.get().set(0.0f);
			}
		}
		
		
		long [][] cropCorr = new long[2][nDim];
		for(i=0;i<nDim;i++)
		{
			cropCorr[0][i]=Math.round((-1)*((temDim[i]-1)*max_fraction_shift));
			cropCorr[1][i]=Math.round((imgDim[i]-1)*max_fraction_shift);
	
		}
		FinalInterval interval = new FinalInterval( cropCorr[0] ,  cropCorr[1] );
		IntervalView< FloatType > ivCCswapped = Views.interval( Views.extendPeriodic( invF1F2 ),interval);
		ImageJFunctions.show(ivCCswapped).setTitle( "General cross-correlation" );
		//ImageJFunctions.show(Views.interval( Views.extendPeriodic( invF1F1I2 ),interval)).setTitle( "Denominator" );
		
		
		//now find max value
		Point shift = new Point(nDim);
		MiscUtils.computeMaxLocation(ivCCswapped,shift);
		for(i=0;i<nDim;i++)
		{
			IJ.log("dim "+Integer.toString(i)+": "+Integer.toString(shift.getIntPosition(i)));
		}
		
		/** final shift of images **/
		long [] shift_pos = new long [nDim];
		//final shift
		shift.localize(shift_pos);
		
		ImagePlus registeredT = ImageJFunctions.wrap(Views.interval(Views.translate(Views.extendZero(template), shift_pos),image),"GNCC registered_template");
		if(nDim>2)
		{
			registeredT.setDimensions(1, (int)ivCCswapped.dimension(2), 1);
		}
		registeredT.show();
			
		
		
	}
	
	public static void calcTermNum(final Img< FloatType > im1, final Img< FloatType > im2,final Img< FloatType > im3, final Img< FloatType > im4)
	{
		//final Img< FloatType > numFin =  im1.copy();
		final Cursor< FloatType > cursFl1 = im1.cursor();
		final Cursor< FloatType > cursFl2 = im2.cursor();
		final Cursor< FloatType > cursFl3 = im3.cursor();
		final Cursor< FloatType > cursFl4 = im4.cursor();
		float temp;
		while(cursFl1.hasNext())
		{
			cursFl1.fwd();
			cursFl2.fwd();
			cursFl3.fwd();
			cursFl4.fwd();
			//t4 = Math.max(cursFl4.get().get(),0);
			//if(Math.abs(t4)>Float.MIN_VALUE)
			{

				temp=cursFl1.get().get()-(cursFl2.get().get()*cursFl3.get().get()/cursFl4.get().get());
				cursFl1.get().set(temp);

			}
		}
		
		
		//return numFin;
	}
	public static void calcTermDenom(final Img< FloatType > im1, final Img< FloatType > im2,final Img< FloatType > im3, final Img< FloatType > im4)
	{
		//final Img< FloatType > numFin =  im1.copy();
		final Cursor< FloatType > cursFl1 = im1.cursor();
		final Cursor< FloatType > cursFl2 = im2.cursor();
		final Cursor< FloatType > cursFl3 = im3.cursor();
		final Cursor< FloatType > cursFl4 = im4.cursor();
		float temp;
		while(cursFl1.hasNext())
		{
			cursFl1.fwd();
			cursFl2.fwd();
			cursFl3.fwd();
			cursFl4.fwd();
			//t4 = Math.max(cursFl4.get().get(),0);
			//if(Math.abs(t4)>Float.MIN_VALUE)
			{

				temp=cursFl1.get().get()-(cursFl2.get().get()*cursFl3.get().get()/cursFl4.get().get());
				cursFl1.get().set(Math.max(temp, 0));

			}
		}
		
		
		//return numFin;
	}
	public static Img< ComplexFloatType > multCompl(final Img< ComplexFloatType > im1, final Img< ComplexFloatType > im2)
	{
		Img< ComplexFloatType > output = im1.copy();
		final Cursor< ComplexFloatType > cursOut = output.cursor();
		final Cursor< ComplexFloatType > cursCompl1 = im1.cursor();
		final Cursor< ComplexFloatType > cursCompl2 = im2.cursor();
		final ComplexFloatType cTemp = new ComplexFloatType();
	    while(cursOut.hasNext())
	    {
	    	cursOut.fwd();
	    	cursCompl1.fwd();
	    	cursCompl2.fwd();
	    	cTemp.set(cursCompl1.get());
	    	cTemp.mul(cursCompl2.get());
	    	cursOut.get().set(cTemp);
	    }
	    
	    return output;
	}
}
