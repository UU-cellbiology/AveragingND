package averagingND;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.IJ;
import net.imglib2.Cursor;
import net.imglib2.FinalDimensions;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
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
import net.imglib2.loops.LoopBuilder;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class MaskedNormCC {
	/** dimensionality of images **/	
	private int nDim;
	
	/** maximum of norm cross-correlation coefficient **/	
	public double dMaxCC;
	
	/** final shift of the template to be aligned (corresponding to CC max) **/
	public long [] dShift;
	
	/** whether to log results to IJ.log window **/
	public boolean bVerbose = true;
	
	/** whether to exclude/ignore voxels that are zeros **/	
	public boolean bZeroMask = false;
	
	/** window limit of CC max localization as a fraction of max displacement **/
	public double [] lim_fractions = null;
	
	/** absolute window limit (in voxels) of CC max localization **/
	public FinalInterval limInterval = null;
	
	/** if false, the limits window of displacement applied to the origin of coordinates,
	 * if true, it is applied around position where both images are centered **/
	public boolean bCenteredLimit = false;
	
	/** the fraction of max displacement where norm CC looks ok,
	 * the boundaries are always strange **/
	final double max_fraction_shift = 0.9;	
	
	ExecutorService es;
	TaskExecutor taskExecutor;
	
	/** 
	 * @param image
	 * @param template
	 */
	public boolean caclulateMaskedNormCC(final RandomAccessibleInterval< FloatType > image, final RandomAccessibleInterval< FloatType > template,  final boolean bShowCC)
	{
		int i;
				
		if(image.numDimensions()!=template.numDimensions())
		{
			IJ.log("Error! Different number of dimensions between reference and template images!");
			return false;
		}
		
		nDim = image.numDimensions();
		if(lim_fractions!= null)
		{
			if(lim_fractions.length!=nDim)
			{
				IJ.log("Error! Constrain on registration for fraction of image has different dimensions!");
				return false;
			}
			if(limInterval != null)
			{
				IJ.log("Warning! multiple constrains (both fraction and pixel)!");
			}

		}
		if(limInterval != null)
		{
			if(limInterval.numDimensions()!=nDim)
			{
				IJ.log("Error! Constrain on registration in pixels has different dimensions!");
				return false;
			}
		}
		
		dShift = new long [nDim];
		long [] imgDim = new long[nDim];
		long [] temDim = new long[nDim];
		image.dimensions(imgDim);
		template.dimensions(temDim);
		final long [] finDim = new long [nDim];
		for(i=0;i<nDim;i++)
		{
			finDim[i]=imgDim[i]+temDim[i]-1;
		}
		
		
		long[] paddedDimensions = new long[nDim];
		long[] fftSize = new long[nDim];
		FFTMethods.dimensionsRealToComplexFast(new FinalDimensions(finDim), paddedDimensions, fftSize);
		
		//padded intervals
		FinalInterval imgIntPad = (FinalInterval) FFTMethods.paddingIntervalCentered(image, new FinalDimensions(paddedDimensions));
		FinalInterval temIntPad = (FinalInterval) FFTMethods.paddingIntervalCentered(template, new FinalDimensions(paddedDimensions));
		
		//padded versions
		IntervalView< FloatType > padImg = Views.interval(Views.extendZero(image),imgIntPad);
		IntervalView< FloatType > padTem = Views.interval(Views.extendZero(template),temIntPad);
		//ImageJFunctions.show(padImg).setTitle( "padded1" );
		//ImageJFunctions.show(padTem).setTitle( "padded2" );
		
		//unity images
		final ArrayImg<FloatType, FloatArray> unityImg = ArrayImgs.floats(imgDim);
		final ArrayImg<FloatType, FloatArray> unityTem = ArrayImgs.floats(temDim);

		if(!bZeroMask)
		{
			unityImg.forEach(t->t.set(1.0f));
			unityTem.forEach(t->t.set(1.0f));
		}
		else
		{
			//reference image
			Cursor< FloatType > unC = unityImg.cursor();
			IntervalView< FloatType > intImg = Views.interval(image,image);
			Cursor< FloatType > inC = intImg.cursor();
			while(unC.hasNext())
			{
				unC.fwd();
				inC.fwd();
				if(inC.get().get()==0.0f)
				{
					unC.get().set(0.0f);
				}
				else
				{
					unC.get().set(1.0f);
				}
			}
			//template image
			unC = unityTem.cursor();
			intImg = Views.interval(template,template);
			inC = intImg.cursor();
			while(unC.hasNext())
			{
				unC.fwd();
				inC.fwd();
				if(inC.get().get()==0.0f)
				{
					unC.get().set(0.0f);
				}
				else
				{
					unC.get().set(1.0f);
				}
			}
		}
		//padded unity images
		IntervalView< FloatType > padUnitImg = Views.interval(Views.extendZero(unityImg),imgIntPad);
		IntervalView< FloatType > padUnitTem = Views.interval(Views.extendZero(unityTem),temIntPad);
		//ImageJFunctions.show(padUnitImg).setTitle( "unipadded1" );
		//ImageJFunctions.show(padUnitTem).setTitle( "unipadded2" );
	
		//squared image values
		//Img< FloatType > sqImg = ImgView.wrap(image).copy();
		//Img< FloatType > sqTem = ImgView.wrap(template).copy();
		//sqImg.forEach(t-> t.mul(t.get()));
		//sqTem.forEach(t-> t.mul(t.get()));
		RandomAccessibleInterval< FloatType > sqImg = Converters.convert(image, ( in,out )-> {out.set(in.get()*in.get());}, new FloatType()); 
		RandomAccessibleInterval< FloatType > sqTem = Converters.convert(template, ( in,out )-> {out.set(in.get()*in.get());}, new FloatType()); 

		//padded squared images
		IntervalView< FloatType > padSqImg = Views.interval(Views.extendZero(sqImg),imgIntPad);
		IntervalView< FloatType > padSqTem = Views.interval(Views.extendZero(sqTem),temIntPad);
		
		//ImageJFunctions.show(padSqImg).setTitle( "sqpadded1" );
		//ImageJFunctions.show(padSqTem).setTitle( "sqpadded2" );
		
		final ImgFactory< ComplexFloatType > factoryComplex = new ArrayImgFactory< ComplexFloatType >(new ComplexFloatType());
		final ImgFactory< FloatType > factoryFloat = new ArrayImgFactory< FloatType >(new FloatType());
		
	
		final int nThreads = Runtime.getRuntime().availableProcessors();
		es = Executors.newFixedThreadPool( nThreads );
		//start with FFT
		final Img< ComplexFloatType > imageFFT2    =   FFT.realToComplex(padImg, factoryComplex,es);
		final Img< ComplexFloatType > templateFFT2 =   FFT.realToComplex(padTem, factoryComplex,es);
		final Img< ComplexFloatType > unImageFFT2  =   FFT.realToComplex(padUnitImg, factoryComplex,es);
		final Img< ComplexFloatType > unTemplateFFT2 = FFT.realToComplex(padUnitTem, factoryComplex,es);
		final Img< ComplexFloatType > sqImageFFT2 =    FFT.realToComplex(padSqImg, factoryComplex,es);
		final Img< ComplexFloatType > sqTemplateFFT2 = FFT.realToComplex(padSqTem, factoryComplex,es);
		
		//conjugates
		FFTMethods.complexConjugate(templateFFT2);	
		FFTMethods.complexConjugate(unTemplateFFT2);	
		FFTMethods.complexConjugate(sqTemplateFFT2);	
		taskExecutor = Parallelization.getTaskExecutor();
		//System.out.println( taskExecutor.suggestNumberOfTasks());
			
		//multiplications
		//reuse already allocated arrays
	    final Img< ComplexFloatType > I1F2F2 = multComplInPlaceSecond(unImageFFT2,sqTemplateFFT2,taskExecutor);
	    final Img< ComplexFloatType > F1F1I2 = multComplInPlaceSecond(unTemplateFFT2,sqImageFFT2,taskExecutor);
	    
		final Img< ComplexFloatType > deNOM = multCompl(unImageFFT2,unTemplateFFT2,taskExecutor);
	    final Img< ComplexFloatType > I1F2 = multComplInPlaceSecond(templateFFT2,unImageFFT2,taskExecutor);
	    final Img< ComplexFloatType > F1I2 = multComplInPlaceSecond(imageFFT2,unTemplateFFT2,taskExecutor);
		final Img< ComplexFloatType > F1F2 = multComplInPlaceSecond(templateFFT2,imageFFT2,taskExecutor);
	        
	    //inverse FFT
		final Img< FloatType > invdeNOM = FFT.complexToReal(deNOM, factoryFloat, new FloatType(),es);	
		final Img< FloatType > invF1F2 = FFT.complexToReal(F1F2, factoryFloat, new FloatType(),es);	
		final Img< FloatType > invF1I2 = FFT.complexToReal(F1I2, factoryFloat, new FloatType(),es);	
		final Img< FloatType > invI1F2 = FFT.complexToReal(I1F2, factoryFloat, new FloatType(),es);	
		final Img< FloatType > invF1F1I2 = FFT.complexToReal(F1F1I2, factoryFloat, new FloatType(),es);	
		final Img< FloatType > invI1F2F2 = FFT.complexToReal(I1F2F2, factoryFloat, new FloatType(),es);
		es.shutdown();
		//ImageJFunctions.show(invdeNOM).setTitle( "denom" );
		//ImageJFunctions.show(invF1F2).setTitle( "invF1F2" );
		//ImageJFunctions.show(invF1I2).setTitle( "invF1I2" );
		//ImageJFunctions.show(invI1F2).setTitle( "invI1F2" );
		//ImageJFunctions.show(invF1F1I2).setTitle( "invF1F1I2" );
		//ImageJFunctions.show(invI1F2F2).setTitle( "invI1F2F2" );
		
		calcTermNum(invF1F2,invF1I2,invI1F2,invdeNOM,taskExecutor);
		calcTermDenom(invF1F1I2,invF1I2,invF1I2,invdeNOM,taskExecutor);
		calcTermDenom(invI1F2F2,invI1F2,invI1F2,invdeNOM,taskExecutor);
		
		
		taskExecutor.close();
		//ImageJFunctions.show(invF1F2).setTitle( "term0" );
		//ImageJFunctions.show(invF1F1I2).setTitle( "term1" );
		//ImageJFunctions.show(invI1F2F2).setTitle( "term2" );

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
			//check if the value within precision tolerance
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
		long [][] cropFraction = new long[2][nDim];

		//determine the maximum area (intervalCrop) of shifts for template with respect to original image
		double nHalfSpan, nCenter;
		
		for(i=0;i<nDim;i++)
		{
			nHalfSpan = 0.5*((temDim[i]-1)+(imgDim[i]-1));
			nCenter = nHalfSpan-(temDim[i]-1);
			cropCorr[0][i] = Math.round(nCenter-nHalfSpan*max_fraction_shift);
			cropCorr[1][i] = Math.round(nCenter+nHalfSpan*max_fraction_shift);	

			if(lim_fractions != null)
			{
				cropFraction[0][i]=Math.round(-nHalfSpan*lim_fractions[i]);
				cropFraction[1][i]=Math.round(nHalfSpan*lim_fractions[i]);				
			}

		}
		
		FinalInterval intervalFull = new FinalInterval(cropCorr[0] ,  cropCorr[1]);
		FinalInterval constrainFr = null;
		if(lim_fractions != null)
		{
			constrainFr = new FinalInterval(cropFraction[0] ,  cropFraction[1]);
		}		

		//Now we need to account for the padding. Since it is changing the origin
		//of coordinates of template with respect to the original image
		long [] nCCOrigin = new long [nDim];
		for(i=0;i<nDim;i++)
		{
			nCCOrigin[i] = (long) ( Math.ceil(0.5*(imgDim[i]))- Math.ceil(0.5*(temDim[i]))); 
		}
		
		
		//In addition, here we do so called FFT "swap quadrants" procedure (or FFT shift)
		//to show frequencies centered/radial. So we can have negative shifts, etc
		//Instead of truly swapping the quadrants we just periodically mirror them (extend periodic)
		IntervalView< FloatType > ivCCswapped;
		//full available shift space
		ivCCswapped =  Views.interval(Views.translate(Views.extendPeriodic( invF1F2 ), nCCOrigin),intervalFull);
		
		//let's apply constrains
			
		//voxel constrains
		if(limInterval != null)
		{
			FinalInterval intConstrainPx; 
			//centered, let's move to the center
			if(bCenteredLimit)
			{
				intConstrainPx = Intervals.intersect(ivCCswapped, Intervals.translate(limInterval, nCCOrigin));
				
			}
			//from the origin, no need to move
			else
			{
				intConstrainPx = Intervals.intersect(ivCCswapped, limInterval);
			}
			ivCCswapped = Views.interval(ivCCswapped, intConstrainPx);
		}
		//fractional constrains
		if(lim_fractions != null)
		{
				//centered, let's move to the center
			if(bCenteredLimit)
			{
				constrainFr =  Intervals.translate(constrainFr, nCCOrigin);
				
			}

			ivCCswapped = Views.interval(ivCCswapped, constrainFr);
		}
		
		
		if(bShowCC)
		{
			ImageJFunctions.show(ivCCswapped).setTitle( "General cross-correlation" );
		}		
		
		//now find max value
		Point shift = new Point(nDim);
		
		FloatType fCCvalue;
	
		fCCvalue = MiscUtils.computeMaxLocation(ivCCswapped,shift);
		
		if (bVerbose)
		{
			String sOutput = "Translation for template to overlap with image (px):\n("+Integer.toString(shift.getIntPosition(0))  ;
			for(i=1;i<nDim;i++)
			{
				sOutput = sOutput +", " +Integer.toString(shift.getIntPosition(i));
				//IJ.log("dim "+Integer.toString(i)+": "+Integer.toString(shift.getIntPosition(i)));
			}
			
			sOutput=sOutput+")\nMaximum normalized cross-correlation value: "+Float.toString(fCCvalue.get());
			IJ.log(sOutput);
		}
		
		// final shift of images 
		shift.localize(dShift);

		// max of cross-correlation
		dMaxCC = fCCvalue.get();

		return true;
	}
	
	public static void calcTermNum(final Img< FloatType > im1, final Img< FloatType > im2,final Img< FloatType > im3, final Img< FloatType > im4, final TaskExecutor te)
	{
		
		LoopBuilder.setImages(im1, im2, im3, im4).multiThreaded(te)
		.forEachPixel(
				(p1,p2,p3,p4)->
				{
					p1.set(p1.getRealFloat()-(p2.getRealFloat()*p3.getRealFloat()/p4.getRealFloat()));
				}
				
				);		

	}
	public static void calcTermDenom(final Img< FloatType > im1, final Img< FloatType > im2,final Img< FloatType > im3, final Img< FloatType > im4, final TaskExecutor te)
	{
		
		LoopBuilder.setImages(im1, im2, im3, im4).multiThreaded(te)
		.forEachPixel(
				(p1,p2,p3,p4)->
				{
					p1.set(Math.max(p1.getRealFloat()-(p2.getRealFloat()*p3.getRealFloat()/p4.getRealFloat()),0));
				}
				
				);	


	}
	public static Img< ComplexFloatType > multCompl(final Img< ComplexFloatType > im1, final Img< ComplexFloatType > im2, final TaskExecutor te)
	{
		Img< ComplexFloatType > output = im1.copy();
		LoopBuilder.setImages(im2,output).multiThreaded(te).forEachPixel((s,t)->t.mul(s));
	    
	    return output;
	}
	public static Img< ComplexFloatType > multComplInPlaceSecond(final Img< ComplexFloatType > im1, final Img< ComplexFloatType > im2, final TaskExecutor te)
	{
		

		LoopBuilder.setImages(im1,im2).multiThreaded(te).forEachPixel((s,t)->t.mul(s));
	    
	    return im2;
	}
}
