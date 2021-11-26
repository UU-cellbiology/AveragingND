package fiji.plugin.RegisterNDFFT;

import ij.IJ;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFT;
import net.imglib2.algorithm.fft2.FFTMethods;
import net.imglib2.algorithm.integral.IntegralImgDouble;
import net.imglib2.converter.RealDoubleConverter;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.FloatArray;
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.iterator.LocalizingIntervalIterator;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class RegisterTranslation
{
	/** dimensionality of images **/	
	int nDim;

	/** number of terms in integral image (summed table) = 2^nDim**/
	int nTerms;
	/** signes of terms in integral image (summed table)**/
	double [] signs;
	/** indexes of min max box coordinates for integrated intensity overlap calculation **/
	int [][] normCCindexes;

	/** final shift of images **/
	long [] shift_pos;
	
	public RegisterTranslation(int nDim)
	{
		this.nDim=nDim;
		
	}
	
	public long[] registerTranslation(final Img< FloatType > image_in, final Img< FloatType > template_in) //throws ImgIOException, IncompatibleTypeException
	{
		
		int i;
		
		double maxShift = 0.5;
		//copies
		final Img< FloatType > image = MiscUtils.copyImage( image_in );
		final Img< FloatType > template = MiscUtils.copyImage( template_in );
		nDim = image.numDimensions();
		shift_pos=new long[nDim];
		long [] imgDim = new long[nDim];
		long [] temDim = new long[nDim];

		image.dimensions(imgDim);
		template.dimensions(temDim);

		// display image and template
		//ImageJFunctions.show( image_in ).setTitle( "input" );
		//ImageJFunctions.show( template_in ).setTitle( "template" );

		
		long nPixImg = imgDim[0];
		long nPixTem = temDim[0];
		
		for(i=1;i<nDim;i++)
		{
			nPixImg *= imgDim[i];
			nPixTem *= temDim[i];			
		}
		
		
		//subtract average values
		MiscUtils.subtractMean(image,nPixImg);
		MiscUtils.subtractMean(template,nPixTem);
		
		
		if(image.numDimensions()!=template.numDimensions())
		{
			IJ.log("different dimensions of input and template!");
			return shift_pos;
		}

		long [][] cropCorr = new long[2][nDim];
		long [][] finDim = new long[2][nDim];
		for(i=0;i<nDim;i++)
		{
			cropCorr[0][i]=Math.round((-maxShift)*temDim[i]);
			cropCorr[1][i]=Math.round(maxShift*imgDim[i]);
			finDim[0][i]=0;
			finDim[1][i]=Math.max(imgDim[i]-1, temDim[i]-1);
		}
				
		long[] paddedDimensions = new long[template.numDimensions()];
		long[] fftSize = new long[template.numDimensions()];
		FFTMethods.dimensionsRealToComplexFast(template, paddedDimensions, fftSize);
		
		final ImgFactory< ComplexFloatType > factoryComplex = new ArrayImgFactory< ComplexFloatType >(new ComplexFloatType());
		final ImgFactory< FloatType > factoryFloat = new ArrayImgFactory< FloatType >(new FloatType());
		
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
		ImageJFunctions.show(invertedT).setTitle( "FFTinv" );
		//swap quadrants (fftshift)
		long [] fftshift = new long [nDim];
		invertedT.max(fftshift);
		//IntervalView< FloatType > ivCC = Views.interval( Views.translate(Views.extendPeriodic( invertedT ),fftshift),interval);
		IntervalView< FloatType > ivCCswapped = Views.interval( Views.extendPeriodic( invertedT ),interval);
		//ImageJFunctions.show(ivCCswapped).setTitle( "nonnorm" );
		
		//calculate summed-area tables for padded images (intensity squared)
		MiscUtils.squareValues(padIm);
		MiscUtils.squareValues(padTem);
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
		
		// let's take a look at one of them
		//ImageJFunctions.show(imgIntegral).setTitle( "integral_image" );
		
	
		final long [] curr_shift = new long [nDim];
		final long [][] overlapImg = new long [2][nDim];
		final long [][] overlapTem = new long [2][nDim];
		
		long [] dimnorm = new long [nDim];
		for (i=0;i<nDim;i++)
		{
			dimnorm[i]=cropCorr[1][i]-cropCorr[0][i]+1;
		}
		

		//calculate indexes and signes for integral image calculation
		calc_index_sign_norm_CC(nDim);

		
		//final correlation normalized image
		ArrayImg<FloatType, FloatArray> normCC = ArrayImgs.floats(dimnorm);
		IntervalView< FloatType > ivCC = Views.translate(normCC, cropCorr[0]);
		final Cursor< FloatType> ivCCRA = ivCC.cursor();
		
		double dIntImg, dIntTem, dFin,dCorr;
		Cursor <FloatType> ivInvCursor = ivCCswapped.cursor();
		while (ivInvCursor.hasNext())
		{
			ivCCRA.fwd();
			ivInvCursor.fwd();
			ivInvCursor.localize(curr_shift);
			calc_overlap(finDim, curr_shift,overlapImg,overlapTem);
			dIntImg=getIntegralIntensityOverlap(imgIntegral,overlapImg);
			dIntTem=getIntegralIntensityOverlap(temIntegral,overlapTem);

			dFin=1.0/Math.sqrt(dIntImg*dIntTem);			
			
			dCorr = ivInvCursor.get().get();
			ivCCRA.get().set((float)(dCorr*dFin));
			
		}
				
		//now find max value
		Point shift = new Point(nDim);
		MiscUtils.computeMaxLocation(ivCC,shift);
		for(i=0;i<nDim;i++)
		{
			IJ.log("dim "+Integer.toString(i)+": "+Integer.toString(shift.getIntPosition(i)));
		}
		//Views.interval(Views.extendZero(template),interval);
		
		shift.localize(shift_pos);
		//Views.interval(Views.translate(Views.extendZero(template), shift_pos),template);
		
		ImageJFunctions.show(ivCC).setTitle( "cross-correlation" );
		ImageJFunctions.show(Views.interval(Views.translate(Views.extendZero(template_in), shift_pos),image_in)).setTitle("registered template");
		return shift_pos;
	}
	
	


	/** function calculates overlap (boxes) of two images of the same size [dims]
	 *  where one of them [overlapTem] is shifted by [shift] with respect to [overlapImg]
	 *   @param dims - dimensions of image (assumed same as template)
	 *   @param shift - shift vector of template
	 *   @param overlapImg - overlap area/volume of image
	 *   @param overlapTem - overlap area/volume of template  
	 * **/
	public static void calc_overlap(long [][] dims, long[] shift, final long [][] overlapImg, final long [][] overlapTem )
	{
		
		for (int i=0;i<shift.length;i++)
		{
			overlapImg[0][i]=Math.max(0, shift[i]);
			overlapImg[1][i]=Math.min(dims[1][i]+1, dims[1][i]+shift[i]+1);
			overlapTem[0][i]=Math.max(0, (-1)*shift[i]);
			overlapTem[1][i]=Math.min(dims[1][i]+1, dims[1][i]-shift[i]+1);

		}
	}
	/** depending on image dimensions,
	 * calculates signs and indexes used to retrieve
	 *  overlap integrated intensity from integral (summed table) image **/
	private void calc_index_sign_norm_CC(int nDim)
	{
		nTerms = (int)Math.pow(2, nDim);
		signs = new double[nTerms];				
		normCCindexes=new int [nTerms][nDim];
		
		int n, nInd, nSum;
		
		//String outS;
		for(byte i=0;i<nTerms;i++)
		{
			//outS="";
			nSum =0 ;
			for(n=0;n<nDim;n++)
			{
				nInd=(i >> (nDim-n-1)) & 1;
				//outS=outS+Integer.toString(nInd);
				normCCindexes[i][n]=nInd;
				nSum+=nInd;
			}
			signs[i]=Math.pow(-1, nDim-nSum);
			//outS=Double.toString(signs[i])+"  "+outS;
			//System.out.println(outS);
		}
		
	}
	
	/** given an integral table and box interval,
	 * calculates integrated intensity over that area **/
	public double getIntegralIntensityOverlap(final RandomAccessibleInterval<DoubleType> integralRAI, final long [][] overlapBox)
	{
		double res = 0;
		int i,j;
		long [] nPos = new long [nDim];
		int nInd =0;
		
		RandomAccess<DoubleType> ra = integralRAI.randomAccess();
		
		for(i=0;i<nTerms; i++)
		{
			for(j=0;j<nDim;j++)
			{
				nInd=normCCindexes[i][j];
				//since integral image has one more pixel along each dimension
				nPos[j]=overlapBox[nInd][j];
			}
			ra.setPosition(nPos);
			res+=ra.get().get()*signs[i];
			
		}
		return res;		
	}
	
	public double honestSum(IntervalView< FloatType > pad, long [][] range)
	{
		double sum = 0;
		long [][] rangex = new long [2][nDim];
		for (int i = 0; i<nDim;i++)
		{
			rangex[0][i]=range[0][i];
			rangex[1][i]=range[1][i]-1;
		}
		IntervalView< FloatType > padRange = Views.interval(pad,rangex[0],rangex[1]);
		for ( final FloatType type : padRange )
		{
			sum += type.get();
			//type.setReal( val*val );
		}
		
		return sum;
		
	}

}
