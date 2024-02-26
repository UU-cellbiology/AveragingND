package registerNDFFT;

import ij.IJ;
import net.imglib2.FinalDimensions;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFT;
import net.imglib2.algorithm.fft2.FFTMethods;
import net.imglib2.converter.ComplexPowerFloatConverter;
import net.imglib2.converter.ComplexPowerGLogFloatConverter;
import net.imglib2.converter.Converters;
import net.imglib2.img.Img;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.LanczosInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.realtransform.PolarToCartesianTransform2D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Scale2D;
import net.imglib2.realtransform.Scale3D;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class RotationCC {
	
	/** dimensionality of images **/	
	private int nDim;
	
	/** 
	 * @param image
	 * @param template
	 */
	public boolean caclulateRotationFFTCC(final RandomAccessibleInterval< FloatType > image, final RandomAccessibleInterval< FloatType > template)
	{
		int i;
		
		//double [] finData;
		

		
		if(image.numDimensions()!=template.numDimensions())
		{
			IJ.log("different dimensions of input and template!");
			return false;
		}
		
		nDim = image.numDimensions();
	
		long [] imgDim = new long[nDim];
		long [] temDim = new long[nDim];
		image.dimensions(imgDim);
		template.dimensions(temDim);
		/*final long [] finDim = new long [nDim];
		for(i=0;i<nDim;i++)
		{
			finDim[i]=imgDim[i]+temDim[i]-1;
		}
		*/
		
		long[] paddedDimensionsIm = new long[nDim];
		long[] paddedDimensionsTem = new long[nDim];
		long[] fftSizeIm = new long[nDim];
		long[] fftSizeTem = new long[nDim];
		FFTMethods.dimensionsRealToComplexFast(new FinalDimensions(imgDim), paddedDimensionsIm, fftSizeIm);
		FFTMethods.dimensionsRealToComplexFast(new FinalDimensions(temDim), paddedDimensionsTem, fftSizeTem);
		
		//System.out.println("done");
		
		//padded intervals
		FinalInterval imgIntPad = (FinalInterval) FFTMethods.paddingIntervalCentered(image, new FinalDimensions(paddedDimensionsIm));
		FinalInterval temIntPad = (FinalInterval) FFTMethods.paddingIntervalCentered(template, new FinalDimensions(paddedDimensionsTem));
		
		//padded versions
		//IntervalView< FloatType > padImg = Views.interval(Views.extendZero(image),imgIntPad);
		//IntervalView< FloatType > padTem = Views.interval(Views.extendZero(template),temIntPad);
		IntervalView< FloatType > padImg = Views.interval(Views.extendMirrorSingle(image),imgIntPad);
		IntervalView< FloatType > padTem = Views.interval(Views.extendMirrorSingle(template),temIntPad);
		//IntervalView< FloatType > padImg = Views.interval(Views.extendBorder(image),imgIntPad);
		//IntervalView< FloatType > padTem = Views.interval(Views.extendBorder(template),temIntPad);
		
		//factories
		final ImgFactory< ComplexFloatType > factoryComplex = new ArrayImgFactory< ComplexFloatType >(new ComplexFloatType());
		final ImgFactory< FloatType > factoryFloat = new ArrayImgFactory< FloatType >(new FloatType());
		
		//start with FFT
		final Img< ComplexFloatType > imageFFT2 = FFT.realToComplex(padImg, factoryComplex);
		final Img< ComplexFloatType > templateFFT2=FFT.realToComplex(padTem, factoryComplex);
		
		ComplexPowerGLogFloatConverter< ComplexFloatType > powerSConv = new ComplexPowerGLogFloatConverter< ComplexFloatType >();
		//ComplexPowerFloatConverter< ComplexFloatType > powerSConv = new ComplexPowerFloatConverter< ComplexFloatType >();
		//Img <FloatType> imagePS =  factoryFloat.create(imageFFT2);
		//ImageJFunctions.show(templateFFT2).setTitle("tempFFTX");
		RandomAccessibleInterval<FloatType> imagePSLog = Converters.convert(( RandomAccessibleInterval< ComplexFloatType > )imageFFT2, powerSConv, new FloatType());
		RandomAccessibleInterval<FloatType> templatePSLog = Converters.convert(( RandomAccessibleInterval< ComplexFloatType > )templateFFT2, powerSConv, new FloatType());


		//ImageJFunctions.show(getLogPolarViewImage(image)).setTitle("image LogPolar");
		ImageJFunctions.show(FFTSwap(imagePSLog)).setTitle("imgFFT");	    
	    ImageJFunctions.show(FFTSwap(templatePSLog)).setTitle("tempFFT");

	    IntervalView<FloatType> transIm = getPolarView(imagePSLog);
	    IntervalView<FloatType> transTem = getPolarView(templatePSLog);

	    //IntervalView<FloatType> transIm = getLogPolarView(imagePSLog);
	    //IntervalView<FloatType> transTem = getLogPolarView(templatePSLog);
	    
		ImageJFunctions.show(transIm).setTitle("imgFFTPolar");	    
	    ImageJFunctions.show(transTem).setTitle("tempFFTPolar");
	    GenNormCC normCC = new GenNormCC();
	    normCC.bExcludeZeros=true;
	    normCC.bZeroX = true;
	    normCC.caclulateGenNormCC(Views.zeroMin(transIm), Views.zeroMin(transTem), 0.25 , false);
	    //normCC.dShift[1]=normCC.dShift[1]-180;
	    System.out.println(Long.toString(normCC.dShift[0]));
	    System.out.println(Long.toString(normCC.dShift[1]));
	    AffineTransform2D rotateTr= new AffineTransform2D();
	   // rotateTr.rotate((normCC.dShift[1]+180)*Math.PI/180.0);
	    rotateTr.rotate((normCC.dShift[1])*Math.PI/180.0);
	    long [] centerPoint = getCenterPoint(template); 
	    long [] centerPointInv =invCenterPoint (centerPoint);
	    IntervalView <FloatType> finrot = Views.interval(
	    									Views.translate(
	    									RealViews.transform(
	    								  		Views.interpolate(
	    								  				Views.translate(
	    								  						Views.extendZero(template),centerPoint),
	    								  			new NLinearInterpolatorFactory<FloatType>()),
	    								  		rotateTr),
	    									centerPointInv),
	    									template);
	    
	    
	    ImageJFunctions.show(finrot).setTitle("registered rotation");
	    return true;
	}
	
	public static IntervalView<FloatType> getPolarView(RandomAccessibleInterval<FloatType> imageIn)
	{
		long [] minR = imageIn.minAsLongArray();
		long [] maxR = imageIn.maxAsLongArray();
		long halfY  = Math.round(0.5*(maxR[1]-1)); 
		maxR[1]=halfY;
		minR[1]=-halfY;
		final FinalInterval fftSQRange = new FinalInterval( minR ,  maxR);	
		IntervalView <FloatType> imageFFTSQ = Views.interval(Views.extendPeriodic(imageIn), fftSQRange);
		PolarToCartesianTransform2D transform = new PolarToCartesianTransform2D();
		//LogPolarToCartesianTransform2D transform = new LogPolarToCartesianTransform2D();
		
		//long radius =(long) Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)));
		long radius = Math.max(maxR[0], maxR[1]);
		
		//Scale2D scaleTransform = new Scale2D(100, Math.PI/180);
		Scale2D scaleTransform = new Scale2D(1.0, Math.PI/180);
		
		FinalInterval polarRange = new FinalInterval( new long []{Math.round(radius*0.25),-89} ,  new long []{radius,90} );
		//Views.extendPeriodic(imagePSLog);
		
	    IntervalView<FloatType> imageFFTRadius= Views.interval(
	    											RealViews.transform(
	    											RealViews.transform(
	    													Views.interpolate(
	    															//Views.translate(Views.extendZero(imageFFTSQ), centerTest), 
	    															Views.extendZero(imageFFTSQ),
	    															//Views.extendBorder(imageFFTSQ),
	    															new LanczosInterpolatorFactory<FloatType>()), 
	    													//new NLinearInterpolatorFactory<FloatType>()), 
	    											transform.inverse()),scaleTransform.inverse()),
	    											polarRange);
	    FinalInterval polarRange2 = new FinalInterval( new long []{Math.round(radius*0.25),-180} ,  new long []{radius,180} );
	    //FinalInterval polarRange2 = new FinalInterval( new long []{0,-180} ,  new long []{1000,180} );
	    return Views.interval(Views.extendPeriodic(imageFFTRadius),polarRange2);
	}
	public static IntervalView<FloatType> getLogPolarView(RandomAccessibleInterval<FloatType> imageIn)
	{
		long [] minR = imageIn.minAsLongArray();
		long [] maxR = imageIn.maxAsLongArray();
		long halfY  = Math.round(0.5*(maxR[1]-1)); 
		maxR[1]=halfY;
		minR[1]=-halfY;
		final FinalInterval fftSQRange = new FinalInterval( minR ,  maxR);	
		IntervalView <FloatType> imageFFTSQ = Views.interval(Views.extendPeriodic(imageIn), fftSQRange);
		//PolarToCartesianTransform2D transform = new PolarToCartesianTransform2D();
		LogPolarToCartesianTransform2D transform = new LogPolarToCartesianTransform2D();
		
		long radius =(long)Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)));
		long radiusLog =(long)Math.ceil(Math.log(Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)))));
		//long radius =(long) Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)));
		double stretchFactor = (double)radiusLog/(double)radius;
		
		//Scale2D scaleTransform = new Scale2D(100, Math.PI/180);
		Scale2D scaleTransform = new Scale2D(stretchFactor, Math.PI/180);
		
		FinalInterval polarRange = new FinalInterval( new long []{0,-89} ,  new long []{radius,90} );
		//Views.extendPeriodic(imagePSLog);
		
	    IntervalView<FloatType> imageFFTRadius= Views.interval(
	    											RealViews.transform(
	    											RealViews.transform(
	    													Views.interpolate(
	    															//Views.translate(Views.extendZero(imageFFTSQ), centerTest), 
	    															Views.extendZero(imageFFTSQ),
	    															//Views.extendBorder(imageFFTSQ),
	    													new NLinearInterpolatorFactory<FloatType>()), 
	    											transform.inverse()),scaleTransform.inverse()),
	    											polarRange);
	    FinalInterval polarRange2 = new FinalInterval( new long []{0,-180} ,  new long []{radius,180} );
	    //FinalInterval polarRange2 = new FinalInterval( new long []{0,-180} ,  new long []{1000,180} );
	    return Views.interval(Views.extendPeriodic(imageFFTRadius),polarRange2);
	}
	
	public static IntervalView<FloatType> getLogPolarViewImage(RandomAccessibleInterval<FloatType> imageIn)
	{
		LogPolarToCartesianTransform2D transform = new LogPolarToCartesianTransform2D();
		//PolarToCartesianTransform2D transform = new PolarToCartesianTransform2D();
		long [] minR = imageIn.minAsLongArray();
		long [] maxR = imageIn.maxAsLongArray();
		long halfX  = Math.round(0.5*(maxR[0]-1)); 
		maxR[0]=halfX;
		minR[0]=-halfX;

		long halfY  = Math.round(0.5*(maxR[1]-1)); 
		maxR[1]=halfY;
		minR[1]=-halfY;
		
		long radius =(long)Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)));
		long radiusLog =(long)Math.ceil(Math.log(Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)))));
		//long radius =(long) Math.ceil(Math.sqrt(Math.pow(maxR[0], 2)+Math.pow(maxR[1], 2)));
		double stretchFactor = (double)radiusLog/(double)radius;
		//Scale2D scaleTransform = new Scale2D(100, Math.PI/180);
		//Scale2D scaleTransform = new Scale2D(1/radius, Math.PI/180);
		Scale2D scaleTransform = new Scale2D(stretchFactor, Math.PI/180);
		
		FinalInterval polarRange = new FinalInterval( new long []{0,-180} ,  new long []{radius,180} );
		//Views.extendPeriodic(imagePSLog);
		
	    IntervalView<FloatType> imageFFTRadius= Views.interval(
	    											RealViews.transform(
	    											RealViews.transform(
	    													Views.interpolate(
	    															Views.translate(Views.extendZero(imageIn), minR), 
	    															//Views.extendZero(imageIn),
	    															//Views.extendBorder(imageFFTSQ),
	    													new NLinearInterpolatorFactory<FloatType>()), 
	    											transform.inverse()),scaleTransform.inverse()),
	    											polarRange);
	    //FinalInterval polarRange2 = new FinalInterval( new long []{0,-180} ,  new long []{radius,180} );
	    //FinalInterval polarRange2 = new FinalInterval( new long []{0,-180} ,  new long []{1000,180} );
	    return imageFFTRadius;
	    //return Views.interval(Views.extendPeriodic(imageFFTRadius),polarRange2);
	}
	
	public static IntervalView<FloatType> FFTSwap(RandomAccessibleInterval<FloatType> imageIn)
	{
		long [] minR = imageIn.minAsLongArray();
		long [] maxR = imageIn.maxAsLongArray();
		long halfY  = Math.round(0.5*(maxR[1]-1)); 
		maxR[1]=halfY;
		minR[1]=-halfY;
		final FinalInterval fftSQRange = new FinalInterval( minR ,  maxR);	
		IntervalView <FloatType> imageFFTSQ = Views.interval(Views.extendPeriodic(imageIn), fftSQRange);
		return imageFFTSQ;
	}
	
	
	public static long [] getCenterPoint (Interval in)
	{
		long [] out = new long [in.numDimensions()];
		double [] minV=in.minAsDoubleArray();
		double [] maxV=in.maxAsDoubleArray();
		for (int i=0;i<in.numDimensions(); i++)
			out[i]=(long)Math.round((-1.0)*(0.5*(maxV[i]-minV[i])+minV[i]));
		return out;
	}
	public static long [] invCenterPoint (long [] centerPoint)
	{
		long [] out = new long [centerPoint.length];
		for (int i=0;i<centerPoint.length; i++)
			out[i]=(-1)*centerPoint[i];
		return out;
	}
}
